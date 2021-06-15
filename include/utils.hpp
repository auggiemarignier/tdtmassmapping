#pragma once

#include "wavetree.h"
#include <string>

std::string enum_to_string(wavetree_perturb_t type);

std::string mkfilename(const char *prefix, const char *file);

std::string mkformatstring(const char *fmt, ...);
