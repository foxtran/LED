# GNU
set(GNU_CXX_WARNINGS     -Wall;-Wextra;-Werror=undef;-Wmissing-include-dirs;-Wpointer-arith;-Winit-self;-Wfloat-equal;-Werror=return-type;-Werror=format=2;-Wredundant-decls;-Wmissing-noreturn;-Wshadow;-Wendif-labels;-Wstrict-aliasing;-Wwrite-strings;-Werror=overflow;-Werror=shift-count-overflow;-Wdate-time;-Wshadow;-Wconversion;-Wuninitialized;-Wnull-dereference;-Wattributes;-Warray-bounds;-Wcast-qual;-Wcast-align;-Wlogical-op;-Wunused-const-variable=2;-Warith-conversion;-Wduplicated-branches;-Wstrict-aliasing=3;-Wsuggest-attribute=const;-Wsuggest-attribute=noreturn;-Wmissing-noreturn;-Wsuggest-attribute=malloc;-Walloca-larger-than=4096)

# Classic Intel
set(Intel_CXX_WARNINGS     -wn=5)

# Intel OneAPI LLVM-based
set(IntelLLVM_CXX_WARNINGS     -Wall)

# Apple Clang
set(AppleClang_CXX_WARNINGS -Wall;-Wextra)
