# QIP-JASMIN Makefile
# 基于JASMIN框架的量子内积计算

# 编译器设置
CXX = g++
CXXFLAGS = -std=c++14 -O2 -Wall -Wextra -DNDIM=2

# JASMIN框架路径 (需要根据实际安装路径调整)
JASMIN_ROOT = /work1/luteng/soft/jusmid
JASMIN_INCLUDE = $(JASMIN_ROOT)/include -I$(JASMIN_ROOT)/include/jasmin
JASMIN_LIB = $(JASMIN_ROOT)/lib

# 包含目录
INCLUDES = -I$(JASMIN_INCLUDE) -I.

# 库文件
LIBS = -L$(JASMIN_LIB) -L$(JASMIN_LIB/jasmin)  -ljasmin -lmpi
#LIBS = -L$(JASMIN_LIB) \
#       -ljasmin \
#       -ljolly \
#       -lboost_system \
#       -lboost_filesystem \
#       -lboost_thread \
#       -lmpi \
#       -lgfortran \
#       -Wl,-rpath,$(JASMIN_LIB)
# 源文件
SOURCES = main_QIP.C QIP.C QIPLevelIntegrator.C
OBJECTS = $(SOURCES:.C=.o)

# 目标可执行文件
TARGET = qip_jasmin

# 默认目标
all: $(TARGET)

# 编译规则
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# 对象文件编译
%.o: %.C
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# 清理
clean:
	rm -f $(OBJECTS) $(TARGET)

# 安装
install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

# 运行测试
test: $(TARGET)
	./$(TARGET) qip.input

# 帮助
help:
	@echo "可用目标:"
	@echo "  all      - 编译程序"
	@echo "  clean    - 清理编译文件"
	@echo "  install  - 安装程序"
	@echo "  test     - 运行测试"
	@echo "  help     - 显示此帮助"

.PHONY: all clean install test help
