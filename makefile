

ECHO_MESSAGE = "Mac OS X"

OBJDIR = obj
DISTDIR = dist
SRCDIR = src
LIBSDIR = libs
INCLUDESDIR = src

CXX = clang++ -v -Wc++11-extensions
CXXFLAGS = -I$(INCLUDESDIR)  -I$(LIBSDIR)/lodepng 

CC = clang++ -v
CFLAGS = $(CXXFLAGS)

APP_NAME = GemTracer

# Find all the C and C++ files we want to compile
# Note the single quotes around the * expressions. The shell will incorrectly expand these otherwise, but we want to send the * directly to the find command.
SOURCES := $(shell find $(SRCDIR) -name '*.cpp' -or -name '*.c')
SOURCES += $(shell find $(LIBSDIR) -name '*.cpp' -or -name '*.c')

# Prepends BUILD_DIR and appends .o to every src file
# As an example, ./your_dir/hello.cpp turns into ./build/./your_dir/hello.cpp.o
OBJS := $(SOURCES:%=$(OBJDIR)/%.o)

LIBS = 

$(info $$SOURCES is [${SOURCES}])
#$(info $$OBJS is [${OBJS}])

##---------------------------------------------------------------------
## BUILD RULES
##---------------------------------------------------------------------

$(OBJDIR):
	mkdir -p $(OBJDIR)


$(DISTDIR):
	mkdir -p $(DISTDIR)


.PHONY: all
all: $(DISTDIR)/$(APP_NAME)
	@echo Build complete for $(ECHO_MESSAGE)


.PHONY: clean
clean:
	yes | rm -rf $(OBJDIR)/* $(DISTDIR)/$(APP_NAME)


$(OBJDIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<


$(OBJDIR)/lodepng.o: $(LIBSDIR)/lodepng/lodepng.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
$(DISTDIR)/$(APP_NAME): $(OBJS)
	mkdir -p $(dir $@)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

