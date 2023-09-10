#include <string>
#include <sys/stat.h>


namespace frizzle
{
namespace filesystem
{
/// create a directory including needed parent directories
bool CreateDirectory(const std::string& path);
/// convert relative paths to absolute ones
std::string FullPath(const std::string& relative_path);
/// yield the directory-part of a string (crops AFTER the last '/')
std::string ParentPath(const std::string& path);
}
}
