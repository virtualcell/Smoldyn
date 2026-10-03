/* HybridGrid: a node-centred Cartesian grid carrying named scalar fields, for
 coupling Smoldyn to an external PDE solver (PDE/particle hybrid simulation).

 The geometry follows VCell's CartesianMesh (vcell-fvsolver), so particle
 rates and histograms line up exactly with VCell's finite volume solver:
   - N nodes per axis; spacing dx = L/(N-1), or L when N == 1;
   - node i sits at x0 + i*dx (the domain centre when N == 1);
   - a point maps to the nearest node, i = (int)((x-x0)*(N-1)/L + 0.5);
   - linear index = i + Nx*(j + Ny*k), i.e. x varies fastest. A numpy array
     of shape (Nz, Ny, Nx) in C order has exactly this layout.

 Header-only and independent of OPTION_VCELL, so the Python molecule histogram
 accessor can use it in every build. */

#ifndef HYBRID_GRID_H
#define HYBRID_GRID_H

#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

class HybridGrid
{
  public:
    HybridGrid(int dim, const double* origin, const double* size, const int* num)
      : dim_(dim)
    {
        if (dim < 1 || dim > 3)
            throw std::invalid_argument("HybridGrid: dim must be 1, 2 or 3");
        for (int d = 0; d < 3; d++) {
            bool used = d < dim;
            origin_[d] = used ? origin[d] : 0.0;
            size_[d] = used ? size[d] : 1.0;
            num_[d] = used ? num[d] : 1;
            if (num_[d] < 1 || size_[d] <= 0)
                throw std::invalid_argument("HybridGrid: num must be >= 1 and size > 0");
            scale_[d] = num_[d] > 1 ? size_[d] / (num_[d] - 1) : size_[d];
        }
    }

    int dim() const { return dim_; }
    const int* num() const { return num_; }
    const double* origin() const { return origin_; }
    const double* size() const { return size_; }
    const double* scale() const { return scale_; }
    long numElements() const { return (long)num_[0] * num_[1] * num_[2]; }

    /* Nearest-node index on one axis, clamped to the grid (VCell asserts instead). */
    int axisIndex(int d, double x) const
    {
        if (num_[d] == 1)
            return 0;
        int i = (int)((x - origin_[d]) * (num_[d] - 1) / size_[d] + 0.5);
        return i < 0 ? 0 : (i >= num_[d] ? num_[d] - 1 : i);
    }

    long index(const double* pos) const
    {
        int i = axisIndex(0, pos[0]);
        int j = dim_ > 1 ? axisIndex(1, pos[1]) : 0;
        int k = dim_ > 2 ? axisIndex(2, pos[2]) : 0;
        return i + (long)num_[0] * (j + (long)num_[1] * k);
    }

    /* Node coordinates for a linear index (VCell CartesianMesh::getVolumeWorldCoord). */
    void center(long idx, double* coords) const
    {
        long rest = idx;
        for (int d = 0; d < 3; d++) {
            long m = rest % num_[d];
            rest /= num_[d];
            coords[d] = num_[d] > 1 ? origin_[d] + m * scale_[d] : origin_[d] + 0.5 * size_[d];
        }
    }

    void setField(const std::string& name, const double* data, long n)
    {
        if (n != numElements())
            throw std::invalid_argument("HybridGrid::setField: field '" + name + "' has " +
                                        std::to_string(n) + " values, grid has " +
                                        std::to_string(numElements()));
        fields_[name].assign(data, data + n);
    }

    /* Stable reference to a field's storage (std::map nodes never move), created
     empty if the field has not been set yet. An empty field evaluates to 0. */
    const std::vector<double>* fieldRef(const std::string& name) { return &fields_[name]; }

    const std::vector<double>* field(const std::string& name) const
    {
        auto it = fields_.find(name);
        return it == fields_.end() ? nullptr : &it->second;
    }

    std::vector<std::string> fieldNames() const
    {
        std::vector<std::string> names;
        for (const auto& kv : fields_)
            if (!kv.second.empty())
                names.push_back(kv.first);
        return names;
    }

    /* Field names referenced by rate expressions (recorded by GridValueProviderFactory). */
    void require(const std::string& name) { required_.insert(name); }
    std::vector<std::string> requiredFields() const
    {
        return std::vector<std::string>(required_.begin(), required_.end());
    }

  private:
    int dim_;
    int num_[3];
    double origin_[3], size_[3], scale_[3];
    std::map<std::string, std::vector<double>> fields_;
    std::set<std::string> required_;
};

#endif
