#pragma once

#include <memory>
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/stdout_color_sinks.h"


namespace polyhedralGravity {

    /**
     * Wrapper Class for spdlog logger
     */
    class PolyhedraleGravityLogger {

        /**
         * The actual spdlog::logger
         */
        const std::shared_ptr<spdlog::logger> _logger;

    public:
        /**
         * Constructs a new PolyhedraleGravityLogger. Further, it registers the new logger in  spdlog's registry with
         * the name POLYHEDRAL_GRAVITY_LOGGER.
         */
        PolyhedraleGravityLogger()
                : _logger(spdlog::stdout_color_mt<spdlog::async_factory>("POLYHEDRAL_GRAVITY_LOGGER")) {
            _logger->set_level(spdlog::level::trace);
        }

        std::shared_ptr<spdlog::logger> getLogger() {
            return _logger;
        }
    };

    /**
     * Single Logger of the Polyhedrale Gravity Model in C++
     * Static initialized
     */
    inline PolyhedraleGravityLogger POLYHEDRAL_GRAVITY_LOGGER{};

}