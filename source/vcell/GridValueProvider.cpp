#include "GridValueProvider.h"

#include <cctype>
#include <cstdlib>
#include <stdexcept>

namespace {

std::string strip(const std::string& s)
{
    std::string out;
    for (char c : s)
        if (!std::isspace((unsigned char)c))
            out += c;
    return out;
}

bool isIdentifier(const std::string& s)
{
    if (s.empty() || !(std::isalpha((unsigned char)s[0]) || s[0] == '_'))
        return false;
    for (char c : s)
        if (!(std::isalnum((unsigned char)c) || c == '_'))
            return false;
    return true;
}

} // namespace

GridValueProvider::GridValueProvider(const std::string& rateExp, std::shared_ptr<HybridGrid> grid)
  : rateExp_(rateExp)
  , grid_(grid)
  , constant_(1.0)
{
    std::string expr = strip(rateExp);
    if (expr.empty())
        throw std::invalid_argument("GridValueProvider: empty rate expression");
    size_t start = 0;
    while (start <= expr.size()) {
        size_t stop = expr.find('*', start);
        std::string token = expr.substr(start, stop == std::string::npos ? std::string::npos : stop - start);
        char* end = nullptr;
        double value = std::strtod(token.c_str(), &end);
        if (!token.empty() && end && *end == '\0') {
            constant_ *= value;
        } else if (isIdentifier(token)) {
            fieldNames_.push_back(token);
            fields_.push_back(grid_->fieldRef(token));
        } else {
            throw std::invalid_argument(
              "GridValueProvider: unsupported rate expression '" + rateExp +
              "' (expected a product of numbers and field names, e.g. 'k*B')");
        }
        if (stop == std::string::npos)
            break;
        start = stop + 1;
    }
}

double GridValueProvider::getConstantValue()
{
    // Smoldyn calls this at load time and treats a throw as "position-dependent".
    if (!fieldNames_.empty())
        throw "rate depends on grid fields";
    return constant_;
}

double GridValueProvider::evaluate(double x, double y, double z)
{
    double pos[3] = { x, y, z };
    long idx = grid_->index(pos);
    double value = constant_;
    for (const std::vector<double>* f : fields_)
        value *= f->empty() ? 0.0 : (*f)[idx];
    return value;
}

double GridValueProvider::getValue(double t, double x, double y, double z, rxnptr rxn)
{
    return evaluate(x, y, z);
}

double GridValueProvider::getValue(double t, double x, double y, double z, rxnptr rxn, char* panelName)
{
    return evaluate(x, y, z);
}

double GridValueProvider::getValue(double t,
  double x,
  double y,
  double z,
  surfactionptr actiondetails,
  char* panelName)
{
    return evaluate(x, y, z);
}

ValueProvider* GridValueProviderFactory::createValueProvider(string& rateExp)
{
    GridValueProvider* provider = new GridValueProvider(rateExp, grid_);
    for (const std::string& name : provider->fieldNames())
        grid_->require(name);
    return provider;
}
