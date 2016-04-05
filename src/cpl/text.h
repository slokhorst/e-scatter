/*!
 * @file src/cpl/text.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__TEXT__HEADER_INCLUDED
#define CPL__TEXT__HEADER_INCLUDED

#include <cinttypes>
#include <string>
#include <vector>

namespace cpl {
namespace text {

/*!
 * @return
 */
inline std::string cat();
template<typename Arg, typename...Args>
std::string cat(const Arg&, const Args&...);
/*!
 * @param[in] file
 * @return
 */
std::string file_to_string(const std::string& file);
/*!
 * @param[in] s
 * @param[in] chars
 * @return
 */
std::string strip_string(const std::string& s, const std::string& chars = " \t\n");
/*!
 * @param[in] s
 * @param[in] c
 * @return
 */
std::vector<std::string> split_string(const std::string& s, char c);
/*!
 * @param[in] s
 * @param[in] c
 * @return
 */
std::string merge_string(const std::vector<std::string>& s, char c);
/*!
 * @param[in] s
 * @param[in] t
 * @return
 */
std::string space_to_tab(const std::string& s, uint8_t t = 4);
/*!
 * @param[in] s
 * @param[in] t
 * @return
 */
std::string tab_to_space(const std::string& s, uint8_t t = 4);
/*!
 * @param[in] x
 * @return
 */
std::string float64(double x);
/*!
 * @param[in] s
 * @return
 */
double float64(const std::string& s);
/*!
 * @param[in] x
 * @return
 */
std::string int32(int32_t x);
/*!
 * @param[in] s
 * @return
 */
int32_t int32(const std::string& s);
/*!
 *
 */
enum color_enum {
	DEFAULT,
	BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, LIGHTGRAY
};
/*!
 * @param[in] s
 * @param[in] c
 * @return
 */
std::string strc(const std::string& s, color_enum c = DEFAULT);
/*!
 * @param[in] s
 * @param[in] c
 * @return
 */
std::string strb(const std::string& s, color_enum c = DEFAULT);

}
}

#include <cpl/text.inl>

#endif