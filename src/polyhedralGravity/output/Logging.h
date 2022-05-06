#pragma once

#include <memory>
#include "spdlog/spdlog.h"
#include "spdlog/async.h"
#include "spdlog/sinks/stdout_color_sinks.h"


namespace polyhedralGravity {

    class PolyhedraleGravityLogger;

    /**
     * Wrapper Class for spdlog logger
     */
    class PolyhedraleGravityLogger {

    public:

        const static PolyhedraleGravityLogger DEFAULT_LOGGER;

    private:

        /**
         * The actual spdlog::logger
         */
        const std::shared_ptr<spdlog::logger> _logger;

    public:
        /**
         * Constructs a new PolyhedraleGravityLogger. Further, it registers the new logger in  spdlog's registry with
         * the name POLYHEDRAL_GRAVITY_LOGGER.
         * TODO async_factory causes bug on Windows OS --> Deadlock on finish (no return 0/ exit)
         */
        PolyhedraleGravityLogger()
                : _logger(spdlog::stdout_color_mt<spdlog::synchronous_factory>("POLYHEDRAL_GRAVITY_LOGGER")) {
            _logger->set_level(spdlog::level::trace);
        }

        [[nodiscard]] std::shared_ptr<spdlog::logger> getLogger() const {
            return _logger;
        }
    };

}