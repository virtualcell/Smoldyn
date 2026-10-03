/* GridValueProvider: position-dependent reaction rates read from fields on a
 HybridGrid that an external PDE solver updates between Smoldyn steps.

 Implements Smoldyn's VCell hybrid hooks (OPTION_VCELL): a reaction rate written
 with a trailing ';' in the configuration file, e.g.

     reaction_cmpt cell r1 A -> 0  kf*B;

 is handed to the ValueProviderFactory. Here the expression must be a product of
 numeric constants and field names ('kf*B', '2.5*B*C', 'B'). Names set with
 'define' are substituted textually by Smoldyn's parser before this point, so
 'kf*B' works when kf is defined; any other name is a grid field. The value
 at a point is the product with each field read at the nearest grid node,
 exactly as vcell-fvsolver's VCellValueProvider does for volume reactions.

 Membrane/surface variants evaluate the same nearest-node volume lookup (the
 panel-indexed membrane lookup used by vcell-fvsolver is not implemented). */

#ifndef GRID_VALUE_PROVIDER_H
#define GRID_VALUE_PROVIDER_H

#include <memory>
#include <string>
#include <vector>

#include "HybridGrid.h"
#include <smoldyn.h>

class GridValueProvider : public ValueProvider
{
  public:
    GridValueProvider(const std::string& rateExp, std::shared_ptr<HybridGrid> grid);
    double getConstantValue();
    double getValue(double t, double x, double y, double z, rxnptr rxn);
    double getValue(double t, double x, double y, double z, rxnptr rxn, char* panelName);
    double getValue(double t, double x, double y, double z, surfactionptr actiondetails, char* panelName);

    const std::vector<std::string>& fieldNames() const { return fieldNames_; }
    double constant() const { return constant_; }

  private:
    double evaluate(double x, double y, double z);

    std::string rateExp_;
    std::shared_ptr<HybridGrid> grid_;
    double constant_;
    std::vector<std::string> fieldNames_;
    std::vector<const std::vector<double>*> fields_;
};

class GridValueProviderFactory : public ValueProviderFactory
{
  public:
    explicit GridValueProviderFactory(std::shared_ptr<HybridGrid> grid) : grid_(grid) {}
    ValueProvider* createValueProvider(string& rateExp);

  private:
    std::shared_ptr<HybridGrid> grid_;
};

/* AbstractMesh over a HybridGrid, used by Smoldyn for 0th-order reactions with
 position-dependent rates (one Poisson draw per grid node, VCell-style). */
class GridMesh : public AbstractMesh
{
  public:
    explicit GridMesh(std::shared_ptr<HybridGrid> grid) : grid_(grid) {}
    void getCenterCoordinates(int volIndex, double* coords) { grid_->center(volIndex, coords); }
    void getDeltaXYZ(double* delta)
    {
        for (int d = 0; d < 3; d++)
            delta[d] = grid_->scale()[d];
    }
    void getNumXYZ(int* num)
    {
        for (int d = 0; d < 3; d++)
            num[d] = grid_->num()[d];
    }

  private:
    std::shared_ptr<HybridGrid> grid_;
};

#endif
