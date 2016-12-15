#pragma once

#include <gsl/gsl>
#include <functional>
#include <vector>
#include <string>
#include "global.hh"


namespace command
{
    class Command_: public std::function<int (
        gsl::span<std::string> const &)>
    {
        public:
            std::string description;

            typedef std::function<int (gsl::span<std::string> const &)> base_type;

            using base_type::base_type;

            template <typename Fn>
            Command_(std::string const &description, Fn fn):
                base_type(fn), description(description) {}
    };

    using Command = Global<Command_>;
}
