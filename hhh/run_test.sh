#!/bin/bash
# QIP-JASMIN 测试脚本

echo "QIP-JASMIN 量子内积计算测试"
echo "============================"

# 编译测试程序
echo "编译测试程序..."
g++ -std=c++14 -O2 -o test_qip test_qip.cpp

if [ $? -eq 0 ]; then
    echo "✓ 编译成功"
else
    echo "✗ 编译失败"
    exit 1
fi

# 运行测试
echo "运行量子内积计算测试..."
./test_qip

if [ $? -eq 0 ]; then
    echo "✓ 测试完成"
else
    echo "✗ 测试失败"
    exit 1
fi

# 清理
echo "清理临时文件..."
rm -f test_qip

echo "测试完成！"
