#pragma once

#include <vector>
#include <memory>
#include "DataSource.h"

/**
 * Interface defining methods giving the calculation some input.
 * This includes the points of interest as well as the source of the data.
 */
class ConfigSource {

public:

    virtual ~ConfigSource() = default;

    /**
     * The vector contains the points for which the polyhedral gravity model should
     * be evaluated.
     * @return vector of points
     */
    virtual std::vector<std::array<double, 3>> getPointsOfInterest() = 0;

    /**
     * The DataSource of the given Polyhedron.
     * @return data source (e. g. a file reader)
     */
    virtual std::shared_ptr<DataSource> getDataSource() = 0;

};

