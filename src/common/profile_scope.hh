/*!
 * @file src/common/profile_scope.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__PROFILE_SCOPE__HEADER_INCLUDED
#define eSCATTER__COMMON__PROFILE_SCOPE__HEADER_INCLUDED

#include <functional>

size_t profile_scope(std::function<void()>);

#endif