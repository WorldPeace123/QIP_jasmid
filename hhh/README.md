# QIP-JASMIN: 基于JASMIN框架的量子内积计算

## 项目简介

本项目将原始的QIP量子内积计算代码适配到JASMIN并行计算框架中，实现了量子计算与高性能并行计算的结合。

## 主要特性

- **量子内积计算**: 使用9量子比特的量子电路模拟器进行内积计算
- **并行计算**: 基于JASMIN框架的分布式并行处理
- **可视化输出**: 支持JaVis数据可视化
- **重启动功能**: 支持计算中断后的恢复

## 文件结构

```
qip-jusmin/
├── main_QIP.C              # 主程序入口
├── QIP.h                   # 量子内积计算类头文件
├── QIP.C                   # 量子内积计算类实现
├── QIPLevelIntegrator.h    # 网格层积分器头文件
├── QIPLevelIntegrator.C    # 网格层积分器实现
├── CMakeLists.txt          # CMake构建文件
├── qip.input              # 输入参数文件
├── vectors.txt            # 输入向量数据
└── README.md              # 说明文档
```

## 编译和运行

### 编译

```bash
mkdir build
cd build
cmake ..
make
```

### 运行

```bash
./qip_jasmin qip.input
```

## 输入参数说明

### 主程序参数 (Main)
- `javis_dump_interval`: 可视化输出间隔
- `javis_dump_dirname`: 可视化输出目录
- `restart_interval`: 重启动文件输出间隔
- `restart_write_dirname`: 重启动文件目录

### 量子计算参数 (QIP)
- `input_vector_file`: 输入向量文件路径
- `output_result_file`: 输出结果文件路径

### 几何参数 (CartesianGeometry)
- `domain_boxes`: 计算域网格
- `x_lo`, `x_hi`: 计算域边界
- `use_ratios`: 是否使用比例
- `periodic_dimensions`: 周期性边界条件

## 算法原理

### 量子内积计算

1. **向量预处理**: 将输入向量进行归一化和符号处理
2. **量子电路构建**: 使用层次归一化方法构建量子电路参数
3. **量子态演化**: 通过9量子比特电路模拟量子态演化
4. **测量**: 使用Swap Test测量量子内积

### 量子电路结构

- **单量子比特门**: X, Z, H, RY门
- **多控制门**: 多控制RY门、X门、Z门
- **双量子比特门**: 受控SWAP门
- **Swap Test**: 用于测量内积的量子算法

## 性能优化

1. **并行化**: 利用JASMIN框架的MPI并行能力
2. **内存管理**: 优化的量子态存储和访问
3. **计算优化**: 高效的复数运算和矩阵操作

## 输出结果

程序会输出以下结果：
- 量子内积计算结果
- 经典内积计算结果
- 平均绝对误差
- 处理时间统计

## 依赖项

- JASMIN框架
- MPI (Message Passing Interface)
- CMake 3.10+
- C++14编译器

## 注意事项

1. 确保JASMIN框架已正确安装和配置
2. 输入向量文件格式为16维向量对
3. 输出结果会追加到指定文件中
4. 支持断点续算功能

## 扩展功能

- 支持不同维度的量子内积计算
- 可扩展的量子电路模板
- 多种量子算法实现
- 分布式量子计算支持
