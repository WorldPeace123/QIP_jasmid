//
// 文件名: QIP.h
// 软件包: JASMIN applications
// 版权  : (c) 2024 量子计算应用
// 描述  : 量子内积计算问题的网格片时间积分算法类.
//

#ifndef included_QIP
#define included_QIP

#include "jasmin/StandardComponentPatchStrategy.h"
#include "jasmin/UniRectangularGridGeometry.h"
#include "jasmin/CellVariable.h"
#include "jasmin/JaVisDataWriter.h"
#include <array>
#include <vector>
#include <cmath>
#include <complex>

using namespace JASMIN;

/**
 * @brief 该类从标准构件网格片策略类 algs::StandardComponentPatchStrategy 派生, 
 * 实现量子内积计算的数值计算子程序.
 *
 * 该类需要从输入文件读取如下参数. 
 *     - \b input_vector_file \n
 *       字符串型, 输入向量文件路径.
 *     - \b output_result_file \n
 *       字符串型, 输出结果文件路径.
 */
class QIP
: public algs::StandardComponentPatchStrategy<NDIM>,
  public tbox::Serializable
{
public:

   /*! @brief 构造函数.
    * @param object_name 输入参数, 字符串, 表示对象名称.
    * @param input_db    输入参数, 指针, 指向输入数据库.
    */
   QIP(const string& object_name,
       tbox::Pointer<tbox::Database> input_db);

   /*!
    * @brief 析构函数.
    */
   virtual ~QIP();

   /// @name 重载基类 algs::StandardComponentPatchStrategy<NDIM> 的函数:
   // @{

   /*!
    * @brief 初始化指定的积分构件.
    */
    void initializeComponent(
            algs::IntegratorComponent<NDIM>* intc) const;

   /**
    * @brief 初始化数据片（支持初值构件）.
    */
   void initializePatchData(
      hier::Patch<NDIM>& patch, 
      const double  time, 
      const bool    initial_time,
      const string& intc_name);

   /*!
    * @brief 完成单个网格片上的数值计算（支持数值构件）.
    */
   void computeOnPatch(hier::Patch<NDIM>& patch,
                             const double  time,
                             const double  dt,
                             const bool    initial_time,
                             const string& intc_name);

   /*! 
    * @brief 在单个网格片上计算稳定性时间步长（支持时间步长构件）.
    */
   double getPatchDt(hier::Patch<NDIM>& patch,
                             const double  time,
                             const bool    initial_time,
                             const int     flag_last_dt,
                             const double  last_dt,
                             const string& intc_name);

   /**
    * @brief 根据边界条件, 填充物理边界影像区数据.
    */
   void setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name);

   //@}
   //
   ///@name 重载基类tbox::Serializable的函数
   //@{

   /*!
   * @brief 将数据成员输出到重启动数据库.
   */
   void putToDatabase(tbox::Pointer<tbox::Database> db);
   
   //@}
   
   ///@name 自定义函数
   //@{

   /*!
    * @brief 注册绘图量.
    */
   void registerPlotData(tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer);

   /*!@brief 量子内积计算核心函数 */
   double q_inner(const std::array<double,16>& a_in, const std::array<double,16>& b_in);

   //@}

private:

   /*!@brief 从输入数据库读入数据.  */
   void getFromInput(tbox::Pointer<tbox::Database> db);

   /*!@brief 从重启动数据库读入数据.  */
   void getFromRestart();

   /*!@brief 注册变量和数据片.  */
   void registerModelVariables();

   /*!@brief 量子电路白盒计算 */
   double Quantum_circuit_whitebox(const std::array<double,15>& theta_a, 
                                   const std::array<double,15>& theta_b,
                                   const std::array<int,16>& flag_a, 
                                   const std::array<int,16>& flag_b);

   /*!@brief 对象名.  */
   string d_object_name;

   /*!@brief 量子态数据片索引号 */
   int d_quantum_state_id;

   /*!@brief 结果数据片索引号 */
   int d_result_id;

   /*!@brief 输入向量文件路径 */
   string d_input_vector_file;

   /*!@brief 输出结果文件路径 */
   string d_output_result_file;

   /*!@brief 处理的数据对数量 */
   int d_processed_pairs;

};

#endif
