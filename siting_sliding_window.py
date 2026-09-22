import os
import arcpy
from arcpy.sa import (
    Raster, FocalStatistics, NbrCircle, SetNull, Con, IsNull, CreateConstantRaster
)

# =============================================================================
# Sliding-window siting with iterative biomass-competition deduction (ArcGIS Pro 3.x)
#
# This script reproduces the ranked candidate-site generation used in the study:
#   1) sum remaining biomethane potential in a nominal 25-km circular window;
#   2) select the cell with the largest window-summed potential;
#   3) if several cells share that maximum, select the one with the higher
#      underlying cell-level biomethane potential;
#   4) assign/remove feedstock inside the selected window; then recalculate.
#
# The production resource raster is 0.01 degree (~1 km nominal resolution), and
# the archived production run used a 25-cell circular neighborhood. Boundary
# cells are handled with DATA so that windows at coastlines/national/raster
# boundaries are truncated to valid cells rather than treated as zero-filled.
#
# Outputs:
#   - outputs/focal/topXXXX_sites.gdb/Sites_Top*_Points
#   - outputs/focal/topXXXX_sites.gdb/Sites_Top*_Buffers
#   - optional focal_sum_####.tif files when SAVE_FOCAL_TIFS = True
# =============================================================================

# -------------------------
# User parameters (edit here)
# -------------------------
IN_RASTER_REL = os.path.join("data", "biomethane_map.tif")
OUT_DIR_REL = os.path.join("outputs", "focal")
N_SITES = 1500
RADIUS_CELLS = 25
SAVE_FOCAL_TIFS = False
# -------------------------


def _project_root():
    # In the ArcGIS Pro Python window, __file__ is not defined; fall back to CWD.
    return os.path.abspath(os.path.dirname(__file__)) if "__file__" in globals() else os.getcwd()


def _safe_delete(path):
    if arcpy.Exists(path):
        try:
            arcpy.management.Delete(path)
        except Exception:
            pass


def run():
    project_root = _project_root()
    in_raster = os.path.join(project_root, IN_RASTER_REL)
    out_dir = os.path.join(project_root, OUT_DIR_REL)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True
    arcpy.env.addOutputsToMap = False
    arcpy.env.parallelProcessingFactor = "0%"

    # Align all derived rasters with the input raster.
    arcpy.env.snapRaster = in_raster
    arcpy.env.cellSize = in_raster
    arcpy.env.extent = in_raster

    out_gdb_name = f"top{N_SITES}_sites.gdb"
    out_gdb = os.path.join(out_dir, out_gdb_name)
    _safe_delete(out_gdb)
    arcpy.management.CreateFileGDB(out_dir, out_gdb_name)
    arcpy.env.workspace = out_gdb
    arcpy.env.scratchWorkspace = out_gdb

    base0 = Raster(in_raster)
    cell_size = float(base0.meanCellWidth)
    radius_map_units = RADIUS_CELLS * cell_size
    sr = arcpy.Describe(in_raster).spatialReference

    print(f"Input raster: {in_raster}")
    print(f"CellSize = {cell_size} | RadiusCells = {RADIUS_CELLS} | RadiusMapUnits = {radius_map_units}")
    print(f"SpatialReference type: {sr.type}")

    sites_fc = os.path.join(out_gdb, f"Sites_Top{N_SITES}_Points")
    buffers_fc = os.path.join(out_gdb, f"Sites_Top{N_SITES}_Buffers")

    arcpy.management.CreateFeatureclass(out_gdb, os.path.basename(sites_fc), "POINT", spatial_reference=sr)
    arcpy.management.AddField(sites_fc, "Iter", "LONG")
    arcpy.management.AddField(sites_fc, "MaxSum", "DOUBLE")

    arcpy.management.CreateFeatureclass(out_gdb, os.path.basename(buffers_fc), "POLYGON", spatial_reference=sr)
    arcpy.management.AddField(buffers_fc, "Iter", "LONG")
    arcpy.management.AddField(buffers_fc, "MaxSum", "DOUBLE")

    # Working rasters (toggle A/B).
    work_a = os.path.join(out_gdb, "working_a")
    work_b = os.path.join(out_gdb, "working_b")
    excl_a = os.path.join(out_gdb, "exclude_a")
    excl_b = os.path.join(out_gdb, "exclude_b")

    try:
        arcpy.management.Copy(in_raster, work_a)
    except Exception:
        Raster(in_raster).save(work_a)

    # Initialize exclusion mask (0 everywhere within current processing extent).
    CreateConstantRaster(0, "INTEGER", cell_size).save(excl_a)

    neighborhood = NbrCircle(RADIUS_CELLS, "CELL")

    # Temporary datasets (overwritten each iteration).
    tmp_search = os.path.join(out_gdb, "tmp_focal_search")
    tmp_tie_underlying = os.path.join(out_gdb, "tmp_tie_underlying")
    tmp_maxcells = os.path.join(out_gdb, "tmp_maxcells")
    tmp_mask = os.path.join(out_gdb, "tmp_mask")
    tmp_pts = os.path.join(out_gdb, "tmp_max_pts")
    tmp_onept = os.path.join(out_gdb, "tmp_onept")
    tmp_buf = os.path.join(out_gdb, "tmp_buf")

    work_in, work_out = work_a, work_b
    excl_in, excl_out = excl_a, excl_b

    # Keep the deduction geometry consistent with the 25-cell focal neighborhood.
    if sr.type == "Geographic":
        buffer_dist = f"{radius_map_units} DecimalDegrees"
    else:
        buffer_dist = f"{RADIUS_CELLS * cell_size}"

    for k in range(1, N_SITES + 1):
        wr = Raster(work_in)
        exr = Raster(excl_in)

        # DATA explicitly ignores NoData neighbors, implementing truncated
        # windows at coastlines/national/raster boundaries.
        focal = FocalStatistics(wr, neighborhood, "SUM", "DATA")

        # Candidate centers must themselves be valid raster cells.
        focal = SetNull(IsNull(wr), focal)

        if SAVE_FOCAL_TIFS:
            focal_out = os.path.join(out_dir, f"focal_sum_{k:04d}.tif")
            focal.save(focal_out)

        # Exclude previously selected areas (exclusion mask == 1 -> NoData).
        focal_search = SetNull(exr, focal, "VALUE = 1")

        for p in [
            tmp_search, tmp_tie_underlying, tmp_maxcells, tmp_mask,
            tmp_pts, tmp_onept, tmp_buf
        ]:
            _safe_delete(p)

        focal_search.save(tmp_search)
        arcpy.management.CalculateStatistics(tmp_search)

        try:
            max_val = float(arcpy.management.GetRasterProperties(tmp_search, "MAXIMUM").getOutput(0))
        except Exception:
            print(f"[{k}] No valid MAXIMUM (all NoData). Stop.")
            break

        if max_val <= 0:
            print(f"[{k}] MAXIMUM <= 0 (max={max_val}). Stop.")
            break

        # Study tie-break rule:
        # among cells sharing the maximum window sum, choose the cell with the
        # highest underlying cell-level biomethane potential.
        tie_underlying = SetNull(Raster(tmp_search) != max_val, wr)
        tie_underlying.save(tmp_tie_underlying)
        arcpy.management.CalculateStatistics(tmp_tie_underlying)
        underlying_max = float(
            arcpy.management.GetRasterProperties(tmp_tie_underlying, "MAXIMUM").getOutput(0)
        )

        selected_cells = Con(
            (Raster(tmp_search) == max_val) & (wr == underlying_max),
            1
        )
        selected_cells.save(tmp_maxcells)
        arcpy.conversion.RasterToPoint(tmp_maxcells, tmp_pts, "VALUE")

        # A second exact tie is extremely unlikely. Use a stable XY ordering only
        # to make such a degenerate case reproducible.
        candidates = []
        with arcpy.da.SearchCursor(tmp_pts, ["SHAPE@XY"]) as cur:
            for row in cur:
                candidates.append(row[0])

        if not candidates:
            print(f"[{k}] Failed to find selected max location point. Stop.")
            break

        x, y = min(candidates, key=lambda xy: (xy[1], xy[0]))
        print(
            f"[{k}] MaxSum={max_val} | Underlying={underlying_max} "
            f"at (x={x}, y={y})"
        )

        pt_geom = arcpy.PointGeometry(arcpy.Point(x, y), sr)

        with arcpy.da.InsertCursor(sites_fc, ["SHAPE@", "Iter", "MaxSum"]) as icur:
            icur.insertRow([pt_geom, k, max_val])

        arcpy.management.CopyFeatures([pt_geom], tmp_onept)
        arcpy.analysis.Buffer(tmp_onept, tmp_buf, buffer_dist, dissolve_option="NONE")

        with arcpy.da.SearchCursor(tmp_buf, ["SHAPE@"]) as cur:
            for row in cur:
                with arcpy.da.InsertCursor(buffers_fc, ["SHAPE@", "Iter", "MaxSum"]) as icur:
                    icur.insertRow([row[0], k, max_val])
                break

        # Polygon -> raster -> binary mask; all feedstock in the selected
        # collection window is removed before the next iteration.
        arcpy.conversion.PolygonToRaster(
            in_features=tmp_buf,
            value_field="OBJECTID",
            out_rasterdataset=tmp_mask,
            cell_assignment="CELL_CENTER",
            cellsize=cell_size,
        )
        mask_bin = Con(IsNull(Raster(tmp_mask)), 0, 1)

        updated_work = Con(mask_bin == 1, 0, wr)
        _safe_delete(work_out)
        updated_work.save(work_out)

        updated_excl = Con(mask_bin == 1, 1, exr)
        _safe_delete(excl_out)
        updated_excl.save(excl_out)

        work_in, work_out = work_out, work_in
        excl_in, excl_out = excl_out, excl_in

        if k % 50 == 0:
            try:
                arcpy.management.ClearWorkspaceCache()
            except Exception:
                pass

    print("Done.")
    print("Points FC :", sites_fc)
    print("Buffers FC:", buffers_fc)
    if SAVE_FOCAL_TIFS:
        print(f"Focal TIFFs: {os.path.join(out_dir, 'focal_sum_0001.tif')} ...")

    arcpy.CheckInExtension("Spatial")


if __name__ == "__main__":
    run()
