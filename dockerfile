FROM python:3.11-slim

# Install system libraries required for GeoPandas, Fiona, GDAL, etc.
RUN apt-get update && apt-get install -y \
    gdal-bin \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    proj-data \
    proj-bin \
    libspatialindex-dev \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Set required GDAL environment variables
ENV GDAL_DATA=/usr/share/gdal
ENV PROJ_LIB=/usr/share/proj

# Copy project files
WORKDIR /app
COPY . /app

# Install Python dependencies
RUN pip install --upgrade pip
RUN pip install \
    flask \
    gunicorn \
    simplekml \
    matplotlib \
    pandas \
    shapely \
    pyproj \
    rtree \
    geopandas \
    fiona

# Expose Render port
ENV PORT=10000

# Start Gunicorn
CMD gunicorn app:app --bind 0.0.0.0:$PORT
