//
// 文件名: QIPLevelIntegrator.h
// 软件包: JASMIN applications
// 版权  : (c) 2024 量子计算应用
// 描述  : 量子内积计算问题的网格层时间积分算法类.
//

#ifndef included_QIPLevelIntegrator
#define included_QIPLevelIntegrator

#include "BasePatchHierarchy.h"
#include "BasePatchLevel.h"
#include "JASMIN_config.h"
#include "tbox/Pointer.h"

#include "StandardComponentPatchStrategy.h"
#include "TimeIntegratorLevelStrategy.h"

#include "CopyIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "NumericalIntegratorComponent.h"
#include "SynchronizeIntegratorComponent.h"

#include "QIP.h"

using namespace JASMIN;

/**
 * @brief 该类从网格层时间积分算法策略类 algs::TimeIntegratorLevelStrategy 派生, 
 * 实现基于网格层的量子内积计算求解流程.
 */
class QIPLevelIntegrator
: public tbox::Serializable,
  public algs::TimeIntegratorLevelStrategy<NDIM>
{
public:

   /*! @brief 构造函数.
    * @param object_name 输入参数, 字符串, 表示对象名称.
    * @param qip_model 输入参数, 指针, 指向量子内积计算模型.
    */
   QIPLevelIntegrator(const string& object_name,
                      QIP* qip_model);

   /*!
    * @brief 析构函数.
    */
   virtual ~QIPLevelIntegrator();

   /// @name 必须实现如下纯虚函数.
   //@{

   /**
    * @brief 初始化网格层时间积分算法.
    *
    * 该函数根据实际应用和计算方法, 创建所有积分构件.
    */
   virtual void initializeLevelIntegrator(
       tbox::Pointer<algs::IntegratorComponentManager<NDIM>> manager);

   /**
    * @brief 初始化指定网格层的数据片.
    * 初值构件类 algs::InitializeIntegratorComponent 支撑该函数的实现.
    *
    * @param hierarchy 输入参数, 指针, 指向待初始化网格层所在的网格片层次结构.
    * @param level_number 输入参数, 整型, 表示待初始化的网格层层号.
    * @param init_data_time 输入参数, 双精度浮点型, 表示初始化的时刻.
    * @param can_be_refined 输入参数, 逻辑型, 表示该网格层可被进一步细化.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param old_level 输入参数, 指针, 指向初始化网格层的旧网格层.
    * @param allocate_data 输入参数, 逻辑型, 真值表示初始化为数据片调度内存空间.
    */
   virtual void initializeLevelData(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int level_number,
       const double init_data_time,
       const bool can_be_refined,
       const bool initial_time,
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> old_level,
       const bool allocate_data);

   /**
    * @brief 计算时间步长.
    * 时间步长构件类 algs::DtIntegratorComponent 支撑该函数的实现.
    *
    * @param level 输入参数, 指针, 指向网格层.
    * @param dt_time 输入参数, 双精度浮点型, 表示当前时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示当前时刻为计算的初始时刻.
    * @param flag_last_dt 输入参数, 整型, 表示上个时间步积分返回的状态.
    * @param last_dt 输入参数, 双精度浮点型, 表示上个时间步长.
    * @return 双精度浮点型, 表示网格层的时间步长.
    */
   virtual double getLevelDt(
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> level,
       const double dt_time,
       const bool initial_time,
       const int flag_last_dt,
       const double last_dt);

   /**
    * @brief 推进网格层一个时间步.
    * 数值积分构件类 algs::NumericalIntegratorComponent 支撑该函数的实现.
    *
    * @param level 输入参数, 指针, 指向待积分的网格层.
    * @param current_time 输入参数, 双精度浮点型, 表示时间步的起始时刻.
    * @param predict_dt 输入参数, 双精度浮点型, 表示为该时间步预测的时间步长.
    * @param max_dt 输入参数, 双精度浮点型, 表示时间步允许的最大时间步长.
    * @param min_dt 输入参数, 双精度浮点型, 表示时间步允许的最小时间步长.
    * @param first_step 输入参数, 逻辑型, 真值当前为重构后或时间步序列的第1步.
    * @param hierarchy_step_number 输入参数, 整型, 表示网格片层次结构的积分步数.
    * @param actual_dt 输出参数, 双精度浮点型, 表示时间步实际采用的时间步长.
    * @return 整型, 表示该时间步积分的状态.
    */
   virtual int advanceLevel(
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> level,
       const double current_time,
       const double predict_dt,
       const double max_dt,
       const double min_dt,
       const bool first_step,
       const int hierarchy_step_number,
       double& actual_dt);

   /**
    * @brief 接受时间相关解.
    * 同步构件类 algs::SynchronizeIntegratorComponent 支撑该函数的实现.
    *
    * @param level 输入参数, 指针, 指向网格层.
    * @param new_time 输入参数, 双精度浮点型, 表示新的时刻.
    * @param deallocate_data 输入参数, 逻辑型, 真值表示接收数值解后, 释放新值数据片的内存空间.
    */
   virtual void acceptTimeDependentSolution(
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> level,
       const double new_time,
       const bool deallocate_data);

   //@}

   /// @name 多层网格应用必须实现的函数
   //@{
   
   /**
    * @brief 是否使用细化时间步进.
    * @return 逻辑型, 真值表示使用细化时间步进.
    */
   virtual bool usingRefinedTimestepping();

   /**
    * @brief 同步新网格片层次结构.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param coarsest_level 输入参数, 整型, 表示最粗网格层层号.
    * @param finest_level 输入参数, 整型, 表示最细网格层层号.
    * @param sync_time 输入参数, 双精度浮点型, 表示同步时刻.
    * @param initial_time 输入参数, 逻辑型, 真值表示初始时刻.
    */
   virtual void synchronizeNewHierarchy(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int coarsest_level,
       const int finest_level,
       const double sync_time,
       const bool initial_time);

   /**
    * @brief 同步较粗网格层.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param coarsest_level 输入参数, 整型, 表示最粗网格层层号.
    * @param sync_time 输入参数, 双精度浮点型, 表示同步时刻.
    * @param old_time 输入参数, 双精度浮点型, 表示旧时刻.
    */
   virtual void synchronizeCoarserLevel(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int coarsest_level,
       const double sync_time,
       const double old_time);

   /**
    * @brief 标记需要细化的网格单元.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param level_number 输入参数, 整型, 表示网格层层号.
    * @param regrid_time 输入参数, 双精度浮点型, 表示重新网格化时刻.
    * @param tag_index 输入参数, 整型, 表示标记索引.
    * @param initial_time 输入参数, 逻辑型, 真值表示初始时刻.
    */
   virtual void tagCellsForRefinement(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int level_number,
       const double regrid_time,
       const int tag_index,
       const bool initial_time);

   //@}

   /// @name 其他必需的虚函数
   //@{
   
   /**
    * @brief 是否使用隐式离散化.
    * @return 逻辑型, 真值表示使用隐式离散化.
    */
   virtual bool usingImplicitDiscretization();

   /**
    * @brief 重置数据到推进前状态.
    * @param level 输入参数, 指针, 指向网格层.
    */
   virtual void resetDataToPreadvanceState(
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> level);

   /**
    * @brief 推进后延长网格层.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param coarsest_level 输入参数, 整型, 表示最粗网格层层号.
    * @param sync_time 输入参数, 双精度浮点型, 表示同步时刻.
    */
   virtual void prolongLevelsAfterAdvance(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int coarsest_level,
       const double sync_time);

   /**
    * @brief 推进前同步网格层.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param coarsest_level 输入参数, 整型, 表示最粗网格层层号.
    * @param sync_time 输入参数, 双精度浮点型, 表示同步时刻.
    */
   virtual void synchronizeLevelsBeforeAdvance(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int coarsest_level,
       const double sync_time);

   /**
    * @brief 标记需要重新映射的网格单元.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param level_number 输入参数, 整型, 表示网格层层号.
    * @param regrid_time 输入参数, 双精度浮点型, 表示重新网格化时刻.
    * @param tag_index 输入参数, 整型, 表示标记索引.
    * @param regrid_start_time 输入参数, 双精度浮点型, 表示重新网格化开始时刻.
    */
   virtual void tagCellsForRemapping(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int level_number,
       const double regrid_time,
       const int tag_index,
       const double regrid_start_time);

   /**
    * @brief 重新映射网格片层.
    * @param level 输入参数, 指针, 指向网格层.
    * @param old_level 输入参数, 指针, 指向旧网格层.
    * @param regrid_time 输入参数, 双精度浮点型, 表示重新网格化时刻.
    */
   virtual void remapPatchLevel(
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> level,
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> old_level,
       const double regrid_time);

   /**
    * @brief 重置网格片层次结构配置.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param coarsest_level 输入参数, 整型, 表示最粗网格层层号.
    * @param finest_level 输入参数, 整型, 表示最细网格层层号.
    */
   virtual void resetHierarchyConfiguration(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int coarsest_level,
       const int finest_level);

   /**
    * @brief 重新网格化.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param coarsest_level 输入参数, 整型, 表示最粗网格层层号.
    * @param finest_level 输入参数, 整型, 表示最细网格层层号.
    * @param regrid_time 输入参数, 双精度浮点型, 表示重新网格化时刻.
    * @return 逻辑型, 表示重新网格化是否成功.
    */
   virtual bool regrid(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int coarsest_level,
       const int finest_level,
       const double regrid_time);

   /**
    * @brief 刷新非均匀工作负载.
    * @param level 输入参数, 指针, 指向网格层.
    * @param workload 输入参数, 数组, 表示工作负载.
    * @param coarsest_level 输入参数, 整型, 表示最粗网格层层号.
    * @param finest_level 输入参数, 整型, 表示最细网格层层号.
    */
   virtual void refreshNonUniformWorkload(
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> level,
       const tbox::Array<double>& workload,
       const int coarsest_level,
       const int finest_level);

   /**
    * @brief 同步较粗网格层（隐式）.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param coarsest_level 输入参数, 整型, 表示最粗网格层层号.
    * @param sync_time 输入参数, 双精度浮点型, 表示同步时刻.
    * @param old_time 输入参数, 双精度浮点型, 表示旧时刻.
    * @return 逻辑型, 表示同步是否成功.
    */
   virtual bool synchronizeCoarserLevelImplicit(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int coarsest_level,
       const double sync_time,
       const double old_time);

   /**
    * @brief 最终化网格层数据.
    * @param hierarchy 输入参数, 指针, 指向网格片层次结构.
    * @param level_number 输入参数, 整型, 表示网格层层号.
    * @param init_data_time 输入参数, 双精度浮点型, 表示初始化时刻.
    * @param can_be_refined 输入参数, 逻辑型, 表示可被细化.
    * @param initial_time 输入参数, 逻辑型, 真值表示初始时刻.
    * @param old_level 输入参数, 指针, 指向旧网格层.
    * @param allocate_data 输入参数, 逻辑型, 真值表示分配数据.
    */
   virtual void finalizeLevelData(
       const tbox::Pointer<hier::BasePatchHierarchy<NDIM>> hierarchy,
       const int level_number,
       const double init_data_time,
       const bool can_be_refined,
       const bool initial_time,
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> old_level,
       const bool allocate_data);

   /**
    * @brief 标记块进行BMR.
    * @param level 输入参数, 指针, 指向网格层.
    * @param tag_array 输入参数, 数组, 表示标记数组.
    */
   virtual void tagBlocksForBMR(
       const tbox::Pointer<hier::BasePatchLevel<NDIM>> level,
       const tbox::Array<int>& tag_array);

   //@}

   /// @name 实现Serializable接口
   //@{
   
   /**
    * @brief 将对象数据写入数据库.
    * @param db 输入参数, 指针, 指向数据库.
    */
   virtual void putToDatabase(tbox::Pointer<tbox::Database> db);

   //@}

private:

   /*!@brief 量子内积计算模型指针 */
   QIP* d_qip_model;

};

#endif

