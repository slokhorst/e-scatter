#pragma once

#include <string>
#include <iostream>
#include <sstream>
#include <experimental/optional>
#include <stdexcept>
#include <map>

#include "fancy_string.hh"
#include "long_string.hh"

/*! \brief Parses command-line options.
 */
namespace argparse
{
    using namespace fancy;

    using std::experimental::optional;
    using std::experimental::nullopt;
    using std::experimental::make_optional;

    inline void format_to(std::ostream &out) {}

    template <typename First, typename ...Rest>
    void format_to(std::ostream &out, First a, Rest &&...rest)
    {
        out << a;
        format_to(out, std::forward<Rest>(rest)...);
    }

    /*! \brief Create a string by outputing arguments to a `std::stringstream`.
     */
    template <typename ...Args>
    std::string format(Args &&...args)
    {
        std::ostringstream ss;
        format_to(ss, std::forward<Args>(args)...);
        return ss.str();
    }

    /*! \brief Exception class for `argparse`.
     */
    class Exception: public std::exception
    {
        std::string msg;

        public:
            Exception(Exception const &other): msg(other.msg) {};
            Exception(Exception &&other) = default;

            template <typename ...Args>
            Exception(Args &&...args):
                msg(format(std::forward<Args>(args)...)) {}

            char const *what() const throw () { return msg.c_str(); }

            ~Exception() throw () {}
    };

    template <typename T>
    optional<T> from_string(std::string const &s)
    {
        std::istringstream iss(s);
        T result;
        iss >> result;

        return (iss.fail() ? nullopt : make_optional(result));
    }


    enum ArgumentKind {
        POSITIONAL,
        FLAG,
        OPTION
    };

    class Args;

    /*! \brief Single option.
     */
    class Argument
    {
        friend class Args;

        std::string tag_;
        std::string description_;
        std::string default_value_;

        ArgumentKind kind_;
        bool optional_;

        public:
            Argument() = default;

            Argument(std::string const &tag, std::string const &description,
                     std::string const &default_value, ArgumentKind kind, 
                     bool optional_ = true):
                tag_(tag), description_(description), default_value_(default_value),
                kind_(kind), optional_(optional_) {}

            /*! \brief Flag argument. */
            Argument(std::string const &tag, std::string const &description):
                tag_(tag),
                description_(description),
                default_value_("0"),
                kind_(FLAG),
                optional_(true) {}

            /*! \brief Valued argument. */
            Argument(std::string const &tag, std::string const &description,
                   std::string const &default_value, bool optional_ = true):
                tag_(tag),
                description_(description),
                default_value_(default_value),
                kind_(OPTION),
                optional_(optional_) {}
    };

    inline Argument positional(
        std::string const &tag, std::string const &description)
    {
        return Argument(tag, description, "", POSITIONAL, false);
    }

    inline Argument option(
        std::string const &tag, std::string const &description,
        std::string const &default_value, bool optional_ = true)
    {
        return Argument(tag, description, default_value, optional_);
    }

    inline Argument option(
        std::string const &tag, std::string const &description, 
        bool optional_ = false)
    {
        return Argument(tag, description, "", optional_);
    }

    inline Argument flag(
        std::string const &tag, std::string const &description)
    {
        return Argument(tag, description);
    }

    /*! \brief Argument parser object, a collection of Options.
     */
    class Args
    {
        std::map<std::string, Argument> option_;
        std::map<std::string, std::string> value_;

        public:
            Args(std::initializer_list<Argument> const &options)
            {
                for (Argument const &o : options)
                    option_[o.tag_] = o;
            }

            std::string usage() const
            {
                std::ostringstream ss;

                auto left_side = [] (std::string const &init)
                {
                    auto line_no = std::make_shared<int>(0);
                    int padding = std::max(0, 18 - int(init.length()));
                    return [line_no, &init, &padding] () -> fancy::ptr
                    {
                        switch (*line_no)
                        {
                            case 0: ++(*line_no);
                                return fancy::compose("  ",
                                    console::fg(colour::fluxom_lime), init,
                                    std::string(padding, ' '),
                                    console::fg(colour::dark_gray), "┃ ",
                                    console::reset());

                            default:
                                return fancy::compose(
                                    console::fg(colour::dark_gray), 
                                    "                     ┃ ", console::reset());
                        }
                    };
                };

                for (auto const &o_ : option_)
                {
                    std::string k;
                    Argument o;
                    std::tie(k, o) = o_;

                    ss << LongString(o.description_, 80, left_side(o.tag_));
                       // << std::endl;
                }

                return ss.str();
            }

            /*! \brief Get the value of an argument. The type of the argument
            should be given as a template parameter. */
            template <typename T>
            optional<T> get(std::string const &tag) const
            {
                return from_string<T>(value_.at(tag));
            }

            /*! \brief Get the value of an argument, return a default value
            if the argument is not present. */
            template <typename T>
            T get_fallback(std::string const &tag, T fallback) const
            {
                optional<T> a = get<T>(tag);
                if (a)
                    return *a;
                else
                    return fallback;
            }

            /*! \brief Parse argument that are given in a range type object.
            The argument should have begin and end methods implemented, i.e.
            be a vector<> or array_view<>.*/
            template <typename Rng>
            void parse(Rng const &args)
            {
                auto begin = args.begin(), end = args.end();

                // start with positionals
                for (auto &o_ : option_)
                {
                    std::string k;
                    Argument o;
                    std::tie(k, o) = o_;

                    if (o.kind_ == POSITIONAL)
                    {
                        if (begin == end)
                            throw Exception("Not enough (positional) arguments.");

                        value_[o.tag_] = *begin++;
                    }
                }

                while (begin != end)
                {
                    std::string tag = *begin;
                    ++begin;

                    if (option_.count(tag) == 0)
                        throw Exception("Unknown option: ", tag);

                    switch  (option_[tag].kind_)
                    {
                        case FLAG:      
                            value_[tag] = "1"; 
                            break;
                        case OPTION:    
                            value_[tag] = *begin++;
                            break;
                        case POSITIONAL: 
                            throw Exception("Positional argument name given as keyword?");
                    }
                }

                for (auto const &kv : option_)
                {
                    if (value_.count(kv.first) == 0)
                    {
                        if (!kv.second.optional_)
                            throw Exception("Missing option: ", kv.first);
                        else
                            value_[kv.first] = kv.second.default_value_;
                    }
                }
            }
    };
}
