import os
import tempfile
import zipfile
import csv
from flask import Flask, render_template, request, send_file, redirect, url_for, flash
from lxml import etree
from shapely.geometry import Polygon, Point, box, MultiPolygon
from shapely.ops import transform
import simplekml
from pyproj import Transformer, CRS

# ----------------------------
# Config
# ----------------------------
UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

app = Flask(__name__)
app.secret_key = "secret123"
app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["MAX_CONTENT_LENGTH"] = 20 * 1024 * 1024  # 20 MB

ALLOWED_EXT = {"kml", "kmz"}


# ----------------------------
# Helpers
# ----------------------------
def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXT


def make_kmz(kml_path, kmz_path):
    """Create a KMZ (zip) with doc.kml inside."""
    with zipfile.ZipFile(kmz_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(kml_path, arcname="doc.kml")


# ----------------------------
# Robust KML parser (manual xml)
# ----------------------------
def parse_kml(filepath):
    """
    Parse KML (or KMZ) and return list of shapely Polygons (lon/lat).
    Supports multiple Polygons and newline-separated coordinates.
    """
    # read bytes (support kmz if file extension .kmz)
    ext = filepath.lower().rsplit(".", 1)[-1]
    if ext == "kmz":
        with zipfile.ZipFile(filepath, "r") as z:
            kml_files = [f for f in z.namelist() if f.lower().endswith(".kml")]
            if not kml_files:
                return []
            with z.open(kml_files[0]) as fh:
                xml = fh.read()
    else:
        with open(filepath, "rb") as fh:
            xml = fh.read()

    root = etree.fromstring(xml)
    ns = {"kml": "http://www.opengis.net/kml/2.2"}

    polygons = []
    coords_nodes = root.findall(".//kml:Polygon//kml:coordinates", ns)

    for node in coords_nodes:
        if node.text is None:
            continue
        raw = node.text.strip()
        parts = raw.replace("\n", " ").replace("\t", " ").split()
        ring = []
        for p in parts:
            vals = p.split(",")
            if len(vals) >= 2:
                try:
                    lon = float(vals[0])
                    lat = float(vals[1])
                    ring.append((lon, lat))
                except Exception:
                    continue
        if len(ring) >= 3:
            polygons.append(Polygon(ring))

    return polygons


# ----------------------------
# Grid generation
# ----------------------------
def make_grid_for_polygon(poly_ll, grid_m, id_start, prefix="G"):
    """
    Create grid for polygon.
    poly_ll: shapely polygon in lon/lat (EPSG:4326)
    grid_m: cell size in meters
    Returns: list of tuples (clipped_polygon_lonlat, area_m2, center_lonlat, id_str), next index
    """
    # project to metric CRS (Web Mercator for simplicity)
    crs_from = CRS.from_epsg(4326)
    crs_to = CRS.from_epsg(3857)
    project_to_m = Transformer.from_crs(crs_from, crs_to, always_xy=True).transform
    project_to_ll = Transformer.from_crs(crs_to, crs_from, always_xy=True).transform

    poly_m = transform(project_to_m, poly_ll)
    minx, miny, maxx, maxy = poly_m.bounds

    cells = []
    idx = id_start
    nx = max(1, int((maxx - minx) / grid_m) + 1)
    ny = max(1, int((maxy - miny) / grid_m) + 1)

    for i in range(nx):
        for j in range(ny):
            x0 = minx + i * grid_m
            y0 = miny + j * grid_m
            cell = box(x0, y0, x0 + grid_m, y0 + grid_m)
            if poly_m.intersects(cell):
                inter_m = poly_m.intersection(cell)
                # area in m2 from metric geometry
                # sometimes intersection returns GeometryCollection -> try to compute area sum
                area_m2 = 0.0
                if inter_m.is_empty:
                    continue
                if isinstance(inter_m, (Polygon, MultiPolygon)):
                    area_m2 = float(inter_m.area)
                else:
                    # geometrycollection: sum polygon areas
                    try:
                        area_m2 = sum([float(g.area) for g in inter_m.geoms if isinstance(g, (Polygon, MultiPolygon))])
                    except Exception:
                        area_m2 = float(inter_m.area) if hasattr(inter_m, "area") else 0.0

                # transform clipped polygon to lon/lat for KML
                try:
                    inter_ll = transform(project_to_ll, inter_m)
                except Exception:
                    # fallback: transform cell centroid only
                    inter_ll = transform(project_to_ll, inter_m)

                center_m = cell.centroid
                center_ll = transform(project_to_ll, center_m)

                id_str = f"{prefix}_{idx:04d}"
                cells.append((inter_ll, area_m2, Point(center_ll.x, center_ll.y), id_str))
                idx += 1

    return cells, idx


# ----------------------------
# Routes
# ----------------------------
@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        # read form
        file = request.files.get("kmlfile")
        grid_size_raw = request.form.get("grid_size")
        prefix = request.form.get("prefix") or "G"
        preview = request.form.get("preview") == "on"
        export_kmz = request.form.get("export_kmz") == "on" if "export_kmz" in request.form else False

        # validations
        if not file or file.filename == "":
            flash("Please upload a KML/KMZ file.", "danger")
            return redirect(request.url)
        if not allowed_file(file.filename):
            flash("Invalid file type. Only .kml or .kmz allowed.", "danger")
            return redirect(request.url)
        try:
            grid_size = float(grid_size_raw)
            if grid_size <= 0:
                raise ValueError("Grid size must be > 0")
        except Exception:
            flash("Grid size must be a positive number (meters).", "danger")
            return redirect(request.url)

        # save upload
        tmpdir = tempfile.mkdtemp(dir=app.config["UPLOAD_FOLDER"])
        filepath = os.path.join(tmpdir, file.filename)
        try:
            file.save(filepath)
        except Exception as e:
            flash(f"Failed to save uploaded file: {e}", "danger")
            return redirect(request.url)

        # parse polygons
        polys = parse_kml(filepath)
        if not polys:
            flash("No polygon found inside KML file.", "danger")
            return redirect(request.url)

        # generate grid KML and collect CSV rows
        kml_out = simplekml.Kml()
        cur_idx = 1
        total_cells = 0
        csv_rows = []  # (Grid_ID, area_ha)

        for i, poly in enumerate(polys, start=1):
            cells, cur_idx = make_grid_for_polygon(poly, grid_size, cur_idx, prefix)
            total_cells += len(cells)

            folder = kml_out.newfolder(name=f"Polygon_{i}_grids")
            # we hide polygon names/labels to prevent map label clutter
            for geom_ll, area_m2, center_ll, gid in cells:
                # add polygon without visible name
                coords = list(geom_ll.exterior.coords)
                pol = folder.newpolygon(name="", outerboundaryis=coords)
                pol.style.polystyle.fill = 0
                pol.style.linestyle.width = 1
                pol.style.linestyle.color = simplekml.Color.red

                # no visible point label (icon hidden and name empty)
                pt = folder.newpoint(name="", coords=[(center_ll.x, center_ll.y)])
                pt.style.iconstyle.scale = 0

                # gather CSV info
                area_ha = area_m2 / 10000.0
                csv_rows.append((gid, round(area_ha, 4)))

        # Save KML
        out_kml = os.path.join(tmpdir, "Grid.kml")
        kml_out.document.description = "Grid generated by web UI. Developed by Krishnasureshkumar FG"
        kml_out.save(out_kml)

        # Save KMZ (if requested) and always create it so users can download
        out_kmz = os.path.join(tmpdir, "Grid.kmz")
        try:
            make_kmz(out_kml, out_kmz)
        except Exception:
            # non-fatal
            out_kmz = None

        # Save CSV attribute table
        out_csv = os.path.join(tmpdir, "Grid_Table.csv")
        try:
            with open(out_csv, "w", newline="", encoding="utf-8") as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(["Grid_ID", "Area_ha"])
                for r in csv_rows:
                    writer.writerow(r)
        except Exception:
            out_csv = None

        # Prepare download URLs
        rel_kml = os.path.basename(tmpdir) + "/Grid.kml"
        rel_kmz = os.path.basename(tmpdir) + "/Grid.kmz" if out_kmz else None
        rel_csv = os.path.basename(tmpdir) + "/Grid_Table.csv" if out_csv else None

        download_kml = url_for("download_file", path=rel_kml)
        download_kmz = url_for("download_file", path=rel_kmz) if rel_kmz else None
        download_csv = url_for("download_file", path=rel_csv) if rel_csv else None

        if preview:
            # show preview and download links
            return render_template(
                "index.html",
                preview_url=download_kml,
                download_url=download_kml,
                download_kmz=download_kmz,
                download_csv=download_csv,
                total_cells=total_cells,
                grid_size=grid_size,
                prefix=prefix
            )

        # if not preview, send KML or KMZ (prefer KMZ if requested)
        if export_kmz and rel_kmz:
            return send_file(os.path.join(app.config["UPLOAD_FOLDER"], rel_kmz), as_attachment=True)
        else:
            return send_file(os.path.join(app.config["UPLOAD_FOLDER"], rel_kml), as_attachment=True)

    # GET
    return render_template("index.html")


# ----------------------------
# Download endpoint
# ----------------------------
@app.route("/download")
def download_file():
    path = request.args.get("path")
    if not path:
        flash("Missing download path.", "danger")
        return redirect(url_for("index"))

    full = os.path.join(app.config["UPLOAD_FOLDER"], path)
    if not os.path.exists(full):
        flash("File not found.", "danger")
        return redirect(url_for("index"))

    return send_file(full, as_attachment=True)


# ----------------------------
# Note: Do NOT use app.run(...) on Render; use gunicorn app:app
# ----------------------------
