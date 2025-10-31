//
// 文件名: main_QIP.C
// 软件包: JASMIN application
// 版权  : (c) 2024 量子计算应用
// 描述  : 主控程序(基于JASMIN框架的量子内积计算).
//
using namespace std;
#include "jasmin/CartesianCoordinates.h"
#include "jasmin/UniRectangularGridGeometry.h"
#include "jasmin/PatchHierarchy.h"
#include "jasmin/HierarchyTimeIntegrator.h"
#include "jasmin/JaVisDataWriter.h"


#include "tbox/JASMINManager.h"
#include "jasmin/tbox/InputManager.h"
#include "jasmin/tbox/RestartManager.h"
#include "jasmin/tbox/InputDatabase.h"
#include "QIPLevelIntegrator.h"
#include "QIP.h"
#include "jaumin/MPI.h"
#include "jaumin/JAUMINManager.h"
#include <iostream>
#include <array>
#include <chrono>
#include <iomanip>
using namespace std;
/*!
 *************************************************************************
 *
 * @brief 基于JASMIN框架的量子内积计算主程序.
 *
 * 该函数分以下几个步骤:
 * -# 预处理: 初始化MPI和JASMIN环境, 解析输入文件, 读取主程序控制参数;
 * -# 创建网格片层次结构和积分算法类对象, 主要包括:
 *    -# 笛卡尔坐标系 geom::CartesianCoordinates<NDIM> ,
 *       均匀矩形网格几何 geom::UniRectangularGridGeometry<NDIM>,
 *       网格片层次结构(单块) hier::PatchHierarchy<NDIM>
 *    -# 网格片积分算法 QIP
 *       - 应用级: 提供量子内积计算的数值计算子程序
 *    -# 网格层积分算法 QIPLevelIntegrator
 *       - 应用级: 提供基于网格层的量子内积计算求解流程
 *    -# 网格片层次结构时间积分算法 algs::HierarchyTimeIntegrator<NDIM>
 * -# 初始化网格片层次结构和物理量数据片;
 * -# 循环: 时间步进;
 * -# 后处理: 释放应用类对象, 释放JASMIN和MPI内部资源.
 *
 ************************************************************************
 */
int main(int argc, char* argv[]) {
  /*******************************************************************************
   *                               预  处  理 *
   *******************************************************************************/
  // 初始化MPI和JASMIN环境.
//  tbox::JASMINManager::startup();
//  tbox::JASMINManager::initialize(argc, argv, 3); 

  // 检查是否为测试模式
  bool test_mode = false;
  if (argc > 1 && string(argv[1]) == "--test") {
    test_mode = true;
    // 在测试模式下，直接运行测试代码，不需要完整的JASMIN初始化
    cout<<"------------------NO JASMIN-----------------\n";
    cout << "QIP量子内积计算测试程序" << endl;
    cout << "========================" << endl;
    
    // 测试向量
    array<double,16> a = {9, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    array<double,16> b = {1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    auto start = chrono::high_resolution_clock::now();
    
    // 创建QIP对象进行计算
    // 创建一个简单的数据库对象用于初始化
    tbox::Pointer<tbox::Database> test_db = new tbox::InputDatabase("TestDB");
    test_db->putString("input_vector_file", "vectors.txt");
    test_db->putString("output_result_file", "qip_results.txt");
    
    QIP* qip_test = new QIP("QIP_Test", test_db);
    double quantum_result = qip_test->q_inner(a, b);
    delete qip_test;
    
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    
    // 经典内积计算
    double classic_result = 0;
    for(int i = 0; i < 16; i++) {
      classic_result += a[i] * b[i];
    }
    
    cout << "输入向量A: ";
    for(int i = 0; i < 16; i++) {
      cout << a[i] << " ";
    }
    cout << endl;
    
    cout << "输入向量B: ";
    for(int i = 0; i < 16; i++) {
      cout << b[i] << " ";
    }
    cout << endl;
    
    cout << "量子内积结果: " << fixed << setprecision(6) << quantum_result << endl;
    cout << "经典内积结果: " << fixed << setprecision(6) << classic_result << endl;
    cout << "绝对误差: " << fixed << setprecision(6) << abs(quantum_result - classic_result) << endl;
    cout << "计算时间: " << duration.count() << " 微秒" << endl;
    
    // 验证结果合理性
    if(abs(quantum_result - classic_result) < 1e-6) {
      cout << "✓ 测试通过: 量子计算结果与经典结果一致" << endl;
    } else {
      cout << "✗ 测试失败: 量子计算结果与经典结果不一致" << endl;
    }
    
    return 0;
  }

  // 非测试模式，正常运行JASMIN框架
  tbox::MPI::init(&argc, &argv);
  tbox::JASMINManager::startup();
  cout<<"------------------Use JAUMIN----------------";
  // 从命令行获取输入文件名.
  string input_filename;
  string restart_read_dirname;
  int restore_num = 0;
  bool is_from_restart = false;

  if ((argc != 2) && (argc != 4)) {
    tbox::pout << "USAGE:  " << argv[0] << " <input filename> "
               << "<restart dir> <restore number> [options]\n"
               << "  options:\n"
               << "  none at this time" << endl;
    tbox::MPI::abort();
    return (-1);
  } else {
    input_filename = argv[1];
    if (argc == 4) {
      restart_read_dirname = argv[2];
      restore_num = atoi(argv[3]);

      is_from_restart = true;
    }
  }

  tbox::plog << "input_filename = " << input_filename << endl;
  tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
  tbox::plog << "restore_num = " << restore_num << endl;




/*
tbox::MPI::init(&argc, &argv);
tbox::JASMINManager::startup();

string input_filename;
string restart_read_dirname;
int restore_num = 0;
bool is_from_restart = false;

if (argc < 2 || argc == 3 || argc > 4) {
    if (tbox::MPI::getRank() == 0) {
        tbox::pout << "USAGE:  " << argv[0] << " <input filename>\n"
                   << "   or:  " << argv[0] << " <input filename> "
                   << "<restart dir> <restore number>" << endl;
    }
    tbox::MPI::abort();
    return -1;
}

input_filename = argv[1];

std::ifstream test_file(input_filename);
if (!test_file.good()) {
    if (tbox::MPI::getRank() == 0) {
        tbox::pout << "Error: cannot open input file: " << input_filename << endl;
    }
    tbox::MPI::abort();
    return -1;
}
test_file.close();

if (argc == 4) {
    restart_read_dirname = argv[2];
    
    try {
        restore_num = std::stoi(argv[3]);
    } catch (const std::exception& e) {
        if (tbox::MPI::getRank() == 0) {
            tbox::pout << "Error: restore number must be an integer" << endl;
        }
        tbox::MPI::abort();
        return -1;
    }
    
    is_from_restart = true;
}

if (tbox::MPI::getRank() == 0) {
    tbox::plog << "input_filename = " << input_filename << endl;
    if (is_from_restart) {
        tbox::plog << "restart_read_dirname = " << restart_read_dirname << endl;
        tbox::plog << "restore_num = " << restore_num << endl;
    }
}

*/


  // 解析输入文件的计算参数到输入数据库, 称之为根数据库.
  tbox::Pointer<tbox::Database> input_db = new tbox::InputDatabase("input_db");
  tbox::InputManager::getManager()->parseInputFile(input_filename, input_db);

  // 从根数据库中获得名称为"Main"的子数据库.
  tbox::Pointer<tbox::Database> main_db = input_db->getDatabase("Main");

  // 从"Main"子数据库中获取可视化输出控制参数.
  int javis_dump_interval =
      main_db->getIntegerWithDefault("javis_dump_interval", 0);
  string javis_dump_dirname = main_db->getString("javis_dump_dirname");

  // 设置重启动参数
  int restart_interval = 0;
  if (main_db->keyExists("restart_interval")) {
    restart_interval = main_db->getInteger("restart_interval");
  }

  string restart_write_dirname;
  if (restart_interval > 0) {
    if (main_db->keyExists("restart_write_dirname")) {
      restart_write_dirname = main_db->getString("restart_write_dirname");
    } else {
      TBOX_ERROR("`restart_interval' > 0, but key `restart_write_dirname'"
                 << "not found in input file.");
    }
  }
  const bool write_restart =
      (restart_interval > 0) && !(restart_write_dirname.empty());

  // 如果是断点续算，打开重启动文件
  tbox::RestartManager* restart_manager = tbox::RestartManager::getManager();
  if (is_from_restart) {
    restart_manager->openRestartFile(restart_read_dirname, restore_num,
                                     tbox::MPI::getNodes());
  }

  {
    /*******************************************************************************
     *                     创建网格片层次结构和积分算法类对象 *
     *******************************************************************************/

    // (1) 创建笛卡尔坐标系.
    tbox::Pointer<geom::CoordinateSystem<NDIM> > coords =
        new geom::CartesianCoordinates<NDIM>();

    // (2) 创建均匀矩形网格几何.
    tbox::Pointer<geom::UniRectangularGridGeometry<NDIM> > grid_geometry =
        new geom::UniRectangularGridGeometry<NDIM>(
            "CartesianGeometry", input_db->getDatabase("CartesianGeometry"),
            coords);

    // (3) 创建网格片层次结构.
    tbox::Pointer<hier::PatchHierarchy<NDIM> > patch_hierarchy =
        new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

    // (4) 创建网格片时间积分算法类（应用级:
    // 提供量子内积计算的数值计算子程序）.
    QIP* qip_model =
        new QIP("QIP", input_db->getDatabase("QIP"));

    // 直接调用q_inner函数实现test_qip.cpp的逻辑算法
    cout << "QIP量子内积计算测试程序" << endl;
    cout << "========================" << endl;
    
    // 测试向量
    array<double,16> a = {1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    array<double,16> b = {1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    auto start = chrono::high_resolution_clock::now();
    
    double quantum_result = qip_model->q_inner(a, b);
    
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    
    // 经典内积计算
    double classic_result = 0;
    for(int i = 0; i < 16; i++) {
      classic_result += a[i] * b[i];
    }
    
    cout << "输入向量A: ";
    for(int i = 0; i < 16; i++) {
      cout << a[i] << " ";
    }
    cout << endl;
    
    cout << "输入向量B: ";
    for(int i = 0; i < 16; i++) {
      cout << b[i] << " ";
    }
    cout << endl;
    
    cout << "量子内积结果: " << fixed << setprecision(6) << quantum_result << endl;
    cout << "经典内积结果: " << fixed << setprecision(6) << classic_result << endl;
    cout << "绝对误差: " << fixed << setprecision(6) << abs(quantum_result - classic_result) << endl;
    cout << "计算时间: " << duration.count() << " 微秒" << endl;
    
    // 验证结果合理性
    if(abs(quantum_result - classic_result) < 1e-6) {
      cout << "✓ 测试通过: 量子计算结果与经典结果一致" << endl;
    } else {
      cout << "✗ 测试失败: 量子计算结果与经典结果不一致" << endl;
    }
    
    // 释放对象并退出程序
    if (qip_model) delete qip_model;
    
    // 释放JASMIN和MPI内部资源
    tbox::JASMINManager::shutdown();
    tbox::MPI::finalize();
    
    return 0;

    // (5) 创建网格层时间积分算法类（应用级:
    // 提供基于网格层的量子内积计算求解流程）.
    QIPLevelIntegrator* qip_level_integrator = new QIPLevelIntegrator(
        "QIPLevelIntegrator", qip_model);

    // (6) 创建网格片层次结构时间积分算法类.
    tbox::Pointer<algs::HierarchyTimeIntegrator<NDIM> > time_integrator =
        new algs::HierarchyTimeIntegrator<NDIM>(
            "HierarchyTimeIntegrator",
            input_db->getDatabase("HierarchyTimeIntegrator"), patch_hierarchy,
            qip_level_integrator);

    // (7) 创建JaVis可视化输出器类, 并注册待输出的绘图量.
    tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_data_writer =
        new appu::JaVisDataWriter<NDIM>("QIP_JaVis_Writer",
                                        javis_dump_dirname);
    qip_model->registerPlotData(javis_data_writer);

    /*******************************************************************************
     *                初 始 化 网 格 片 层 次 结 构 和 物 理 量 *
     *******************************************************************************/
    tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
    tbox::pout << "initialize hierarchy ... " << endl;
    tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;

    // 初始化网格片层次结构和物理量.
    time_integrator->initializeHierarchy();

    // 关闭重启动输入文件
    if (is_from_restart) restart_manager->closeRestartFile();

    // 输出初始时刻数据场.
    javis_data_writer->writePlotData(patch_hierarchy,
                                     time_integrator->getIntegratorStep(),
                                     time_integrator->getIntegratorTime());

    /************************************************************************************
     *                              时  间  步  进 *
     ************************************************************************************/

    double loop_time = time_integrator->getIntegratorTime();
    double loop_time_end = time_integrator->getEndTime();
    int iteration_num = time_integrator->getIntegratorStep();

    while ((loop_time < loop_time_end) && time_integrator->stepsRemaining()) {
      tbox::pout << "\n\n++++++++++++++++++++++++++++++++++++++++++++" << endl;
      tbox::pout << "Step # " << iteration_num << endl;
      tbox::pout << "begin time = " << loop_time << endl;

      // 将数值解推进一个时间步, 返回其时间步长.
      double dt_actual = time_integrator->advanceHierarchy();

      loop_time += dt_actual;
      iteration_num++;

      // 输出重启动数据
      if ((write_restart) && (iteration_num % restart_interval == 0)) {
        tbox::RestartManager::getManager()->writeRestartFile(
            restart_write_dirname, iteration_num);
      }

      // 输出可视化数据.
      if ((javis_dump_interval > 0) &&
          (iteration_num % javis_dump_interval) == 0) {
        javis_data_writer->writePlotData(patch_hierarchy, iteration_num,
                                         loop_time);
      }

      tbox::pout << "dt = " << dt_actual << endl;
      tbox::pout << "end time = " << loop_time << endl;
      tbox::pout << "++++++++++++++++++++++++++++++++++++++++++++\n\n" << endl;
    }

    /************************************************************************************
     *                               模  拟  结  束 *
     ************************************************************************************/
    // 释放对象.
    if (qip_level_integrator) delete qip_level_integrator;
    if (qip_model) delete qip_model;
  }

  // 释放JASMIN和MPI内部资源.
  tbox::JASMINManager::shutdown();
  tbox::MPI::finalize();

  return (0);
}
