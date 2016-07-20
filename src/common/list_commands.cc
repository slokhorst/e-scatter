#include "fancy_string.hh"
#include "long_string.hh"
#include "command.hh"

using namespace fancy;
using console::fg;
using console::reset;
using command::Command;

Command list_commands("list",
    "Give a list of the different commands contained in this little program.",
    [] (gsl::span<std::string> const &args) -> int
{
    std::cout << "test ground, list of commands:\n";
    for (auto const &p : Command::dir())
    {
        int line_no = 0;
        auto left_side = [&line_no, &p] () -> fancy::ptr
        {
            switch (line_no)
            {
                case 0: ++line_no;
                    return fancy::compose(
                        fg(colour::fluxom_lime), p.first,
                        fg(colour::dark_gray), ":  ", reset());

                default:
                    return fancy::compose(
                        fg(colour::dark_gray), "  ┃ ", reset());
            }
        };

        std::cout << LongString(p.second->description, 80, left_side) << std::endl;
    }

    return 0;
});
