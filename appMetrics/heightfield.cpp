#include "heightfield.h"

#include <algorithm>
#include <cmath>

#include <vector>
#include <cstdint>
#include <limits>

#include <iostream>

const double HeightField::flat = 1.0e-8;

/*!
\brief Create a flat heightfield.
\param box Rectangle domain of the terrain.
\param nx, ny Samples.
\param v Constant elevation.
*/
HeightField::HeightField(const Box2& box, int nx, int ny, const double& v)
    : ScalarField2(box, nx, ny, v)
{
}

/*!
\brief Create a heightfield.
\param box Rectangle domain of the terrain.
\param nx, ny Samples.
\param v Set of elevation values.
*/
HeightField::HeightField(const Box2& box, int nx, int ny, const std::vector<double>& v)
    : ScalarField2(box, nx, ny, v)
{
}

/*!
\brief Create a heightfield from a scalar field.
This constructor provides implicit conversion.
\param s Scalar field.
*/
HeightField::HeightField(const ScalarField2& s)
    : ScalarField2(s)
{
}

/*!
\brief Create a heightfield from an image.
\param box Rectangle domain of the terrain.
\param image Elevation image.
\param a, b Minimum and maximum elevation range.
\param grayscale Boolean set to false if the image is provided in color.
*/
HeightField::HeightField(const Box2& box, const QImage& image, const double& a, const double& b, bool grayscale)
    : ScalarField2(box, image, a, b, grayscale)
{
}

/*!
\brief Compute the elevation of a given position on the terrain.
\param p Point.
\param triangular Boolean, use triangular interpolation if set to true, bilinear interpolation otherwise.
*/
inline double HeightField::Height(const Vector2& p, bool triangular) const
{
    double u, v;
    int i, j;
    cellCoords(p, i, j, u, v);

    double z = 0.0;
    if (isValidCell(i, j)) {
        if (triangular) {
            if (u > v) {
                z = (1.0 - u) * at(i, j) + (u - v) * at(i + 1, j) + v * at(i + 1, j + 1);
            }
            else {
                z = (1.0 - v) * at(i, j) + u * at(i + 1, j + 1) + (v - u) * at(i, j + 1);
            }
        }
        else {
            z = Math::Bilinear(at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1), u, v);
        }
    }
    return z;
}


/*!
\brief Compute the vertex position on the terrain.
\param p Point.
\param triangular Boolean, use triangular interpolation if set to true, bilinear interpolation otherwise.
*/
Vector3 HeightField::Vertex(const Vector2& p, bool triangular) const
{
    return Vector3(p[0], p[1], Height(p, triangular));
}


/*!
\brief Compute the gradient at a given array vertex.

\param i,j Integer coordinates of the array vertex.
*/
Vector2 HeightField::Gradient(int i, int j) const
{
    Vector2 n;

    // Gradient along x axis
    if (i == 0)
        n[0] = (at(i + 1, j) - at(i, j)) * inverseCellSize[0];
    else if (i == nx - 1)
        n[0] = (at(i, j) - at(i - 1, j)) * inverseCellSize[0];
    else
        n[0] = (at(i + 1, j) - at(i - 1, j)) * 0.5 * inverseCellSize[0];

    // Gradient along y axis
    if (j == 0)
        n[1] = (at(i, j + 1) - at(i, j)) * inverseCellSize[1];
    else if (j == ny - 1)
        n[1] = (at(i, j) - at(i, j - 1)) * inverseCellSize[1];
    else
        n[1] = (at(i, j + 1) - at(i, j - 1)) * 0.5 * inverseCellSize[1];

    return n;
}

Vector2 HeightField::Gradient(const Vector2 &p) const
{
    int i, j;
    double u, v;

    cellCoords(p, i, j, u, v);

    if (!isValidCell(i, j)) return Vector2(0,0);

    Vector2 g00 = Gradient(i, j);
    Vector2 g10, g01, g11;
    if (isValidCell(i+1, j))   g10 = Gradient(i+1, j);
    else                             g10 = g00;
    if (isValidCell(i, j+1))   g01 = Gradient(i, j+1);
    else                             g01 = g00;
    if (isValidCell(i+1, j+1)) g11 = Gradient(i+1, j+1);
    else                             g11 = g00;

    Vector2 g0 = (1-v)*g00 + v*g01;
    Vector2 g1 = (1-v)*g10 + v*g11;
    return (1-u)*g0 + u*g1;
}


/*!
\brief Compute the normal for a given position on the terrain.

Note that this function may be expensive to compute.

\param p Point.
\param triangular Boolean, use triangle normals if set to true, bilinear interpolation of normals at vertices otherwise.
*/
Vector3 HeightField::Normal(const Vector2& p, bool triangular) const
{
    double u, v;
    int i, j;
    cellCoords(p, i, j, u, v);

    // Test position
    if (!isValidCell(i, j))
        return Vector3(0,0,0);

    if (triangular) {
        if (u > v) {
            return Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).Normal();
        }
        else {
            return Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).Normal();
        }
    }
    else {
        return Normalized(Bilinear(Normal(i, j), Normal(i + 1, j), Normal(i + 1, j + 1), Normal(i, j + 1), u, v));
    }
}


/*!
\brief Compute the normal at a given sample.

This function uses the weighted sum (area) of the normals of the
triangles sharing the point on the grid. The returned vector is normalized.

\param i,j Integer coordinates of the sample.
*/
Vector3 HeightField::Normal(int i, int j) const
{
  Vector3 n;
  if (i == 0)
  {
    if (j == 0)
    {
      // Corner: 0/1
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() + Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal();
    }
    else if (j == ny - 1)
    {
      // Corner: 5
      n = Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
    else
    {
      // Edge: 0/1/5
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() + Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal()
        + Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
  }
  else if (i == nx - 1)
  {
    if (j == 0)
    {
      // Corner: 2
      n = Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal();

    }
    else if (j == ny - 1)
    {
      // Corner: 3/4
      n = Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal() + Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal();
    }
    else
    {
      // Edge: 2/3/4
      n = Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal();
    }
  }
  else
  {
    if (j == 0)
    {
      // Edge: 0/1/2
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal();
    }
    else if (j == ny - 1)
    {
      // Edge: 3/4/5
      n = Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
    else
    {
      // Face: 0/1/2/3/4/5
      n = Triangle(Vertex(i, j), Vertex(i + 1, j), Vertex(i + 1, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i + 1, j + 1), Vertex(i, j + 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j + 1), Vertex(i - 1, j)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j), Vertex(i - 1, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i - 1, j - 1), Vertex(i, j - 1)).AreaNormal() +
        Triangle(Vertex(i, j), Vertex(i, j - 1), Vertex(i + 1, j)).AreaNormal();
    }
  }
  return Normalized(n);
}

void HeightField::thermalErode(int iterations, double talusAngleDeg, double rate, bool useDiagonals)
{
    if (nx <= 1 || ny <= 1 || iterations <= 0) return;
    rate = std::clamp(rate, 0.0, 1.0);

    const double pi = 3.14159265358979323846;
    const double tanTalus = std::tan(talusAngleDeg * pi / 180.0);

    const double dx = cellSize[0];
    const double dy = cellSize[1];

    struct Off { int di, dj; double dist; };
    std::vector<Off> neigh;
    neigh.reserve(useDiagonals ? 8 : 4);
    neigh.push_back({ 1, 0, dx });
    neigh.push_back({-1, 0, dx });
    neigh.push_back({ 0, 1, dy });
    neigh.push_back({ 0,-1, dy });
    if (useDiagonals) {
        const double dd = std::hypot(dx, dy);
        neigh.push_back({ 1, 1, dd });
        neigh.push_back({ 1,-1, dd });
        neigh.push_back({-1, 1, dd });
        neigh.push_back({-1,-1, dd });
    }
   std::vector<double> delta(nx * ny, 0.0);
    for (int it = 0; it < iterations; ++it) {
        std::fill(delta.begin(), delta.end(), 0.0);

        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                const int id = cellId(i, j);
                const double h = field[id];

                int nIds[8];
                double excess[8];
                int nCount = 0;

                for (const auto& o : neigh) {
                    const int ni = i + o.di;
                    const int nj = j + o.dj;
                    if (!isValidCell(ni, nj)) continue;

                    const int nid = cellId(ni, nj);
                    const double hn = field[nid];

                    const double threshold = tanTalus * o.dist;
                    const double ex = (h - hn) - threshold;
                    if (ex > 0.0) {
                        nIds[nCount] = nid;
                        excess[nCount] = ex;
                        ++nCount;
                    }
                }

                if (nCount == 0) continue;

                for (int k = 0; k < nCount; ++k) {
                    const double move = 0.5 * rate * excess[k];
                    delta[id]      -= move;
                    delta[nIds[k]] += move;
                }
            }
        }

        for (int id = 0; id < nx * ny; ++id)
            field[id] += delta[id];
    }
}

void HeightField::hydraulicErode(
    int iterations,
    double rain,
    double evaporation,
    double flowRate,
    double capacityK,
    double erosionRate,
    double depositionRate,
    bool useDiagonals
)
{
    if (nx <= 1 || ny <= 1 || iterations <= 0) return;

    rain            = std::max(0.0, rain);
    evaporation     = std::clamp(evaporation, 0.0, 1.0);
    flowRate        = std::max(0.0, flowRate);
    capacityK       = std::max(0.0, capacityK);
    erosionRate     = std::max(0.0, erosionRate);
    depositionRate  = std::max(0.0, depositionRate);

    const int N = nx * ny;
    const double eps = 1.0e-12;

    std::vector<double> water(N, 0.0);
    std::vector<double> sed(N,   0.0);

    std::vector<double> waterNext(N, 0.0);
    std::vector<double> sedNext(N,   0.0);

    struct Off { int di, dj; double dist; };
    std::vector<Off> neigh;
    neigh.reserve(useDiagonals ? 8 : 4);

    const double dx = cellSize[0];
    const double dy = cellSize[1];
    neigh.push_back({ 1, 0, dx });
    neigh.push_back({-1, 0, dx });
    neigh.push_back({ 0, 1, dy });
    neigh.push_back({ 0,-1, dy });

    if (useDiagonals) {
        const double dd = std::hypot(dx, dy);
        neigh.push_back({ 1, 1, dd });
        neigh.push_back({ 1,-1, dd });
        neigh.push_back({-1, 1, dd });
        neigh.push_back({-1,-1, dd });
    }

    for (int it = 0; it < iterations; ++it)
    {
        for (int id = 0; id < N; ++id)
            water[id] += rain;

        std::fill(waterNext.begin(), waterNext.end(), 0.0);
        std::fill(sedNext.begin(),   sedNext.end(),   0.0);

        for (int i = 0; i < nx; ++i)
        for (int j = 0; j < ny; ++j)
        {
            const int id = cellId(i, j);
            const double hw = field[id] + water[id];

            int bestId = id;
            double bestHW = hw;
            double bestDist = 1.0;

            for (const auto& o : neigh) {
                const int ni = i + o.di;
                const int nj = j + o.dj;
                if (!isValidCell(ni, nj)) continue;
                const int nid = cellId(ni, nj);
                const double nhw = field[nid] + water[nid];
                if (nhw < bestHW) {
                    bestHW = nhw;
                    bestId = nid;
                    bestDist = o.dist;
                }
            }

            if (bestId == id) {
                waterNext[id] += water[id];
                sedNext[id]   += sed[id];
                continue;
            }

            const double dhw = hw - bestHW;
            if (dhw <= 0.0) {
                waterNext[id] += water[id];
                sedNext[id]   += sed[id];
                continue;
            }

            const double flow = std::min(water[id], flowRate * dhw);

            const double sedMoved = sed[id] * (flow / (water[id] + eps));

            waterNext[id]     += (water[id] - flow);
            sedNext[id]       += (sed[id]   - sedMoved);

            waterNext[bestId] += flow;
            sedNext[bestId]   += sedMoved;

            const double slope = (dhw / (bestDist + eps));   
            const double capacity = capacityK * flow * slope; 

            if (sedNext[id] > capacity) {
                const double dep = depositionRate * (sedNext[id] - capacity);
                sedNext[id]   -= dep;
                field[id]     += dep;
            } else {
                double ero = erosionRate * (capacity - sedNext[id]);

                const double maxErode = 0.05 * std::min(cellSize[0], cellSize[1]); 
                ero = std::min(ero, maxErode);

                const double floorH = 0.0; 
                ero = std::min(ero, field[id] - floorH);


                field[id]     -= ero;
                sedNext[id]   += ero;
            }
        }

        water.swap(waterNext);
        sed.swap(sedNext);

        const double keep = 1.0 - evaporation;
        for (int id = 0; id < N; ++id)
            water[id] *= keep;
    }
}

static inline double clamp01(double x) { return std::max(0.0, std::min(1.0, x)); }
static inline double deg2rad(double d) { return d * 3.14159265358979323846 / 180.0; }

static inline uint32_t hash_u32(uint32_t x)
{
    x ^= x >> 16;
    x *= 0x7feb352dU;
    x ^= x >> 15;
    x *= 0x846ca68bU;
    x ^= x >> 16;
    return x;
}

ScalarField2 HeightField::slopeField(bool useDiagonals) const
{
    ScalarField2 out(getDomain(), getSizeX(), getSizeY(), 0.0);
    const double dx = cellSize[0];
    const double dy = cellSize[1];
    const double eps = 1e-12;

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    {
        auto H = [&](int x, int y) -> double {
            x = std::max(0, std::min(nx-1, x));
            y = std::max(0, std::min(ny-1, y));
            return field[cellId(x,y)];
        };

        double hx = (H(i+1,j) - H(i-1,j)) / (2.0*dx + eps);
        double hy = (H(i,j+1) - H(i,j-1)) / (2.0*dy + eps);

        if (useDiagonals) {
            double hxy1 = (H(i+1,j+1) - H(i-1,j-1)) / (2.0*std::hypot(dx,dy) + eps);
            double hxy2 = (H(i+1,j-1) - H(i-1,j+1)) / (2.0*std::hypot(dx,dy) + eps);
            hx = 0.7*hx + 0.15*hxy1 + 0.15*hxy2;
            hy = 0.7*hy + 0.15*hxy1 - 0.15*hxy2;
        }

        out(i,j) = std::sqrt(hx*hx + hy*hy);
    }
    return out;
}

ScalarField2 HeightField::drainageAreaD8() const
{
    const double dx = cellSize[0];
    const double dy = cellSize[1];
    struct Off { int di, dj; double dist; };
    const Off neigh[8] = {
        { 1, 0, dx }, {-1, 0, dx }, { 0, 1, dy }, { 0,-1, dy },
        { 1, 1, std::hypot(dx,dy) }, { 1,-1, std::hypot(dx,dy) },
        {-1, 1, std::hypot(dx,dy) }, {-1,-1, std::hypot(dx,dy) }
    };

    const int N = nx*ny;
    std::vector<int> to(N, -1);
    std::vector<int> indeg(N, 0);
    std::vector<double> area(N, 1.0);

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    {
        const int id = cellId(i,j);
        const double h0 = field[id];

        int best = -1;
        double bestDrop = 0.0;

        for (const auto& o : neigh) {
            int ni = i + o.di, nj = j + o.dj;
            if (!isValidCell(ni,nj)) continue;
            int nid = cellId(ni,nj);
            double drop = (h0 - field[nid]) / (o.dist + 1e-12);
            if (drop > bestDrop) { bestDrop = drop; best = nid; }
        }

        to[id] = best; 
        if (best != -1) indeg[best] += 1;
    }

    std::vector<int> q;
    q.reserve(N);
    for (int id = 0; id < N; ++id) if (indeg[id] == 0) q.push_back(id);

    for (size_t qi = 0; qi < q.size(); ++qi) {
        int u = q[qi];
        int v = to[u];
        if (v != -1) {
            area[v] += area[u];
            indeg[v] -= 1;
            if (indeg[v] == 0) q.push_back(v);
        }
    }

    ScalarField2 A(getDomain(), getSizeX(), getSizeY(), 0.0);

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
        A(i,j) = area[cellId(i,j)];

    return A;
}

ScalarField2 HeightField::wetnessIndex(const ScalarField2& A, const ScalarField2& s) const
{

    ScalarField2 w(getDomain(), getSizeX(), getSizeY(), 0.0);

    const double eps = 1e-12;
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    {
        double a = A.at(i,j);
        double sl = s.at(i,j);
        w(i,j) = std::log((a + 1.0) / (sl + eps));
    }
    return w;
}

ScalarField2 HeightField::lightField(double sunAzimuthDeg, double sunAltitudeDeg) const
{
    const double az = deg2rad(sunAzimuthDeg);
    const double al = deg2rad(sunAltitudeDeg);

    const double sx = std::cos(al) * std::cos(az);
    const double sy = std::cos(al) * std::sin(az);
    const double sz = std::sin(al);

    ScalarField2 L(getDomain(), getSizeX(), getSizeY(), 0.0);

    const double dx = cellSize[0];
    const double dy = cellSize[1];
    const double eps = 1e-12;

    auto H = [&](int x, int y) -> double {
        x = std::max(0, std::min(nx-1, x));
        y = std::max(0, std::min(ny-1, y));
        return field[cellId(x,y)];
    };

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    {
        double hx = (H(i+1,j) - H(i-1,j)) / (2.0*dx + eps);
        double hy = (H(i,j+1) - H(i,j-1)) / (2.0*dy + eps);

        double nxv = -hx, nyv = -hy, nzv = 1.0;
        double invLen = 1.0 / std::sqrt(nxv*nxv + nyv*nyv + nzv*nzv + eps);
        nxv *= invLen; nyv *= invLen; nzv *= invLen;

        double dot = nxv*sx + nyv*sy + nzv*sz;
        L(i,j) = clamp01(dot);
    }
    return L;
}

static inline double smoothstep(double e0, double e1, double x)
{
    double t = (x - e0) / (e1 - e0);
    t = std::max(0.0, std::min(1.0, t));
    return t*t*(3.0 - 2.0*t);
}

ScalarField2 HeightField::vegetationSuitability(const ScalarField2& slope,
                                                const ScalarField2& wetness,
                                                const ScalarField2& light,
                                                double slopeMaxDeg,
                                                double wetnessMu,
                                                double wetnessSigma) const
{

    ScalarField2 P(getDomain(), getSizeX(), getSizeY(), 0.0);


    const double slopeMax = std::tan(deg2rad(slopeMaxDeg));
    const double inv2sig2 = 1.0 / (2.0 * wetnessSigma * wetnessSigma + 1e-12);

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    {
        double s = slope.at(i,j);
        double Rs = 1.0 - smoothstep(0.6*slopeMax, slopeMax, s);
        double w = wetness.at(i,j);
        double Rw = std::exp(-(w - wetnessMu)*(w - wetnessMu) * inv2sig2);

        double Rl = clamp01(light.at(i,j));

        double p = std::min(Rs, std::min(Rw, Rl));
        P(i,j) = clamp01(p);
    }
    return P;
}

void HeightField::vegetationSeed(ScalarField2& V, const ScalarField2& P, double seedProb, uint32_t seed) const
{
    seedProb = std::max(0.0, std::min(1.0, seedProb));
    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    {
        uint32_t h = hash_u32(seed ^ uint32_t(cellId(i,j)*2654435761u));
        double u = (h & 0x00FFFFFFu) / double(0x01000000u); // [0,1)
        double p = seedProb * clamp01(P.at(i,j));
        if (u < p) V(i,j) = std::min(1.0, V(i,j) + 0.25);
    }
}

void HeightField::vegetationStep(ScalarField2& V,
                                 const ScalarField2& P,
                                 double growth,
                                 double death,
                                 double diffusion,
                                 double dt) const
{
    growth = std::max(0.0, growth);
    death  = std::max(0.0, death);
    diffusion = std::max(0.0, diffusion);

    ScalarField2 Vn(getDomain(), getSizeX(), getSizeY(), 0.0);


    auto getV = [&](int x, int y) -> double {
        x = std::max(0, std::min(nx-1, x));
        y = std::max(0, std::min(ny-1, y));
        return V(x,y);
    };

    for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j)
    {
        double v = V(i,j);
        double p = clamp01(P.at(i,j));

        double dv = (growth * p) * v * (1.0 - v) - (death * (1.0 - p)) * v;

        double lap = getV(i+1,j) + getV(i-1,j) + getV(i,j+1) + getV(i,j-1) - 4.0*v;

        double vn = v + dt*dv + diffusion * lap;
        Vn(i,j) = clamp01(vn);
    }

    V = Vn;
}

double heuristic(int x1, int y1, int x2, int y2, double dx, double dy) {
    return std::sqrt(std::pow((x1 - x2) * dx, 2) + std::pow((y1 - y2) * dy, 2));
}

double HeightField::computeSlopeCost(int x1, int y1, int x2, int y2, double slopeWeight) const
{
    double d2 = std::sqrt(std::pow((x1 - x2) * cellSize[0], 2) + 
                          std::pow((y1 - y2) * cellSize[1], 2));

    double h1 = at(x1, y1);
    double h2 = at(x2, y2);
    double dh = std::abs(h2 - h1);

    double slope = dh / d2;

    double slopePenalty = 1.0 + (slopeWeight * slope * slope);
    
    if (slope > 1.0) slopePenalty *= 10.0;

    return d2 * slopePenalty;
}

std::vector<Vector3> HeightField::findRoadPath(const Vector3& startWorld, const Vector3& endWorld, double slopeWeight) const
{
    std::vector<Vector3> path;

    Index2 startIdx = cellCoords(Vector2(startWorld));
    Index2 endIdx = cellCoords(Vector2(endWorld));

    int startX = startIdx.x();
    int startY = startIdx.y();
    int endX = endIdx.x();
    int endY = endIdx.y();

    if (!isValidCell(startX, startY) || !isValidCell(endX, endY)) {
        return path;
    }

    int w = nx;
    int h = ny;
    
    std::vector<std::vector<double>> costSoFar(w, std::vector<double>(h, std::numeric_limits<double>::max()));
    std::vector<std::vector<Index2>> cameFrom(w, std::vector<Index2>(h, Index2(-1, -1)));
    
    std::priority_queue<PathNode, std::vector<PathNode>, std::greater<PathNode>> openSet;

    openSet.push({startX, startY, 0.0, heuristic(startX, startY, endX, endY, cellSize[0], cellSize[1]), 0.0, -1, -1});
    costSoFar[startX][startY] = 0.0;

    while (!openSet.empty()) {
        PathNode current = openSet.top();
        openSet.pop();

        if (current.x == endX && current.y == endY) {
            break;
        }

        if (current.gCost > costSoFar[current.x][current.y]) {
            continue;
        }

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                if (dx == 0 && dy == 0) continue;

                int nx = current.x + dx;
                int ny = current.y + dy;

                if (isValidCell(nx, ny)) {
                    double moveCost = computeSlopeCost(current.x, current.y, nx, ny, slopeWeight);
                    double newCost = current.gCost + moveCost;

                    if (newCost < costSoFar[nx][ny]) {
                        costSoFar[nx][ny] = newCost;
                        double h = heuristic(nx, ny, endX, endY, cellSize[0], cellSize[1]);
                        openSet.push({nx, ny, newCost, h, newCost + h, current.x, current.y});
                        cameFrom[nx][ny] = Index2(current.x, current.y);
                    }
                }
            }
        }
    }

    Index2 current = Index2(endX, endY);
    if (cameFrom[endX][endY].x() == -1) {
        return path;
    }

    while (current.x() != -1) {
        path.push_back(Vertex(current.x(), current.y()));
        current = cameFrom[current.x()][current.y()];
        
        if (current.x() == startX && current.y() == startY) {
            path.push_back(Vertex(startX, startY));
            break;
        }
    }

    std::reverse(path.begin(), path.end());
    
    return path;
}

ScalarField2 HeightField::computeSnowCover(double snowAltitude, double steepnessThreshold) const
{
    ScalarField2 snow(getDomain(), nx, ny, 0.0);

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            double h = at(i, j);

            if (h < snowAltitude)
            {
                snow(i, j) = 0.0;
                continue;
            }

            Vector3 n = Normal(i, j);

            double slope = 1.0 - n[2]; 

            if (slope > steepnessThreshold)
            {
                snow(i, j) = 0.0;
                continue;
            }

            double altitudeFactor = (h - snowAltitude) / 200.0;
            altitudeFactor = std::clamp(altitudeFactor, 0.0, 1.0);

            double slopeFactor = 1.0 - (slope / steepnessThreshold);
            slopeFactor = std::clamp(slopeFactor, 0.0, 1.0);

            snow(i, j) = altitudeFactor * slopeFactor;
        }
    }
    return snow;
}

void HeightField::applyBrush(const Vector3& p, double radius, double intensity)
{
    Index2 center = cellCoords(Vector2(p));
    
    int R = static_cast<int>(std::ceil(radius / cellSize[0]));

    for (int i = -R; i <= R; ++i) {
        for (int j = -R; j <= R; ++j) {
            
            int nx = center.x() + i;
            int ny = center.y() + j;

            if (isValidCell(nx, ny)) {
                
                double dist = std::sqrt(i*i + j*j) * cellSize[0];

                if (dist <= radius) {

                    double falloff = 0.5 * (1.0 + std::cos(3.14159265359 * dist / radius));

                    int id = cellId(nx, ny);
                    field[id] += intensity * falloff;
                }
            }
        }
    }
}

static double pointSegmentDistanceSq(const Vector2& p, const Vector2& a, const Vector2& b, double& t) {
    Vector2 ab = b - a;
    Vector2 ap = p - a;
    double lenSq = ab[0]*ab[0] + ab[1]*ab[1];
    if (lenSq < 1e-8) {
        t = 0.0;
        return ap[0]*ap[0] + ap[1]*ap[1];
    }
    t = (ap[0] * ab[0] + ap[1] * ab[1]) / lenSq;
    t = std::clamp(t, 0.0, 1.0);
    Vector2 closest = a + ab * t;
    Vector2 diff = p - closest;
    
    return diff[0]*diff[0] + diff[1]*diff[1];
}

void HeightField::applyRoadTerraforming(const std::vector<Vector3>& path, 
                                        double roadWidth, 
                                        double blendWidth,
                                        ScalarField2* vegetation)
{
    if (path.size() < 2) return;

    double totalRadius = roadWidth + blendWidth;
    
    for (size_t k = 0; k < path.size() - 1; ++k) {
        Vector3 pA_world = path[k];
        Vector3 pB_world = path[k+1];

        Index2 iA = cellCoords(Vector2(pA_world));
        Index2 iB = cellCoords(Vector2(pB_world));
        Vector2 vA(iA.x(), iA.y());
        Vector2 vB(iB.x(), iB.y());

        int minX = std::min(iA.x(), iB.x()) - (int)totalRadius;
        int maxX = std::max(iA.x(), iB.x()) + (int)totalRadius;
        int minY = std::min(iA.y(), iB.y()) - (int)totalRadius;
        int maxY = std::max(iA.y(), iB.y()) + (int)totalRadius;

        minX = std::max(0, minX); maxX = std::min(nx - 1, maxX);
        minY = std::max(0, minY); maxY = std::min(ny - 1, maxY);

        for (int i = minX; i <= maxX; ++i) {
            for (int j = minY; j <= maxY; ++j) {
                Vector2 currentP(i, j);

                double t;
                double distSq = pointSegmentDistanceSq(currentP, vA, vB, t);
                double dist = std::sqrt(distSq);

                double distMeters = dist * cellSize[0];

                if (distMeters < totalRadius) {

                    
                    double targetHeight = pA_world[2] * (1.0 - t) + pB_world[2] * t;                    
                    double alpha = 0.0; 

                    if (distMeters < roadWidth) {

                        alpha = 1.0;
                        
                        if (vegetation) (*vegetation)(i, j) = 0.0; 

                    } else {
                        double relativeDist = (distMeters - roadWidth) / blendWidth;
                        
                        alpha = 0.5 * (1.0 + std::cos(3.14159265359 * relativeDist));
                        
                        if (vegetation) (*vegetation)(i, j) *= (1.0 - alpha);
                    }

                    int id = cellId(i, j);
                    double currentHeight = field[id];
                    field[id] = alpha * targetHeight + (1.0 - alpha) * currentHeight;
                }
            }
        }
    }
}