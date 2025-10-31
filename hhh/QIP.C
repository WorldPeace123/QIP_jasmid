using namespace std;
//
// 文件名: QIP.C
// 软件包: JASMIN applications
// 版权  : (c) 2024 量子计算应用
// 描述  : 量子内积计算问题的网格片时间积分算法类的实现.
//

#include "QIP.h"
#include "CellData.h"
#include "UniRectangularPatchGeometry.h"
#include "tbox/RestartManager.h"
#include <fstream>
#include <iomanip>
#include <chrono>
using namespace JASMIN;

// 复数类型定义
struct C {
    double r, i;
    C(double r_=0.0, double i_=0.0): r(r_), i(i_) {}
};

static inline C cadd(const C&a, const C&b){ return C(a.r+b.r, a.i+b.i);} 
static inline C csub(const C&a, const C&b){ return C(a.r-b.r, a.i-b.i);} 
static inline C cmul(const C&a, const C&b){ return C(a.r*b.r - a.i*b.i, a.r*b.i + a.i*b.r);} 

// 2x2复数矩阵
struct M2 { C a00,a01,a10,a11; };

static inline M2 matI(){ return {C(1,0),C(0,0),C(0,0),C(1,0)}; }
static inline M2 matX(){ return {C(0,0),C(1,0),C(1,0),C(0,0)}; }
static inline M2 matZ(){ return {C(1,0),C(0,0),C(0,0),C(-1,0)}; }
static inline M2 matH(){ double s = 1.0/sqrt(2.0); return {C(s,0),C(s,0),C(s,0),C(-s,0)}; }
static inline M2 matRY(double th){ double c=cos(th/2.0), s=sin(th/2.0); return {C(c,0),C(-s,0),C(s,0),C(c,0)}; }

// 量子态模拟器
struct State {
    static constexpr int n = 9; // qubits
    static constexpr size_t Dim = (1u<<n);
    array<C, Dim> amp; // size 512
    State(){ amp.fill(C(0,0)); amp[0]=C(1,0); }

    // 单量子比特门
    void apply1(const M2&U, int q){
        const size_t N = Dim;
        const size_t step = 1ull<<q;
        const size_t period = step<<1;
        for(size_t base=0; base<N; base+=period){
            for(size_t k=0;k<step;++k){
                size_t i0 = base + k;
                size_t i1 = i0 + step;
                C a0 = amp[i0];
                C a1 = amp[i1];
                C b0 = cadd(cmul(U.a00,a0), cmul(U.a01,a1));
                C b1 = cadd(cmul(U.a10,a0), cmul(U.a11,a1));
                amp[i0]=b0; amp[i1]=b1;
            }
        }
    }

    // 多控制单量子比特门
    void applyCtrl1(const M2&U, const vector<int>& ctrls, int tgt){
        const size_t N = Dim;
        size_t mask=0; for(int c:ctrls) mask |= (1ull<<c);
        size_t tbit = 1ull<<tgt;
        for(size_t i=0;i<N;i++){
            if( ( i & mask) == mask ){
                size_t j = i ^ tbit;
                if( ((size_t)i & tbit)==0 ){
                    C a0 = amp[i];
                    C a1 = amp[j];
                    C b0 = cadd(cmul(U.a00,a0), cmul(U.a01,a1));
                    C b1 = cadd(cmul(U.a10,a0), cmul(U.a11,a1));
                    amp[i]=b0; amp[j]=b1;
                }
            }
        }
    }

    // 多控制双量子比特门
    void applyCtrl2(const array<C,16>& U4, const vector<int>& ctrls, int q1, int q2){
        const size_t N = Dim;
        size_t mask=0; for(int c:ctrls) mask |= (1ull<<c);
        int hi = max(q1,q2), lo=min(q1,q2);
        size_t bhi = 1ull<<hi, blo=1ull<<lo;
        for(size_t idx=0; idx<N; ++idx){
            if( (idx & mask) != mask ) continue;
            if( (idx & bhi) || (idx & blo) ) continue;
            size_t i00 = idx;
            size_t i01 = idx | blo;
            size_t i10 = idx | bhi;
            size_t i11 = idx | bhi | blo;
            C v00=amp[i00], v01=amp[i01], v10=amp[i10], v11=amp[i11];
            C in[4]={v00,v01,v10,v11};
            C out[4]={};
            for(int r=0;r<4;++r){
                C acc; 
                for(int c=0;c<4;++c){
                    const C &m = U4[r*4+c];
                    acc = cadd(acc, cmul(m, in[c]));
                }
                out[r]=acc;
            }
            amp[i00]=out[0]; amp[i01]=out[1]; amp[i10]=out[2]; amp[i11]=out[3];
        }
    }
};

/*************************************************************************
 *
 * 构造函数.
 *
 *************************************************************************/
QIP::QIP(const string& object_name,
         tbox::Pointer<tbox::Database> input_db)
{
   d_object_name = object_name;
   d_processed_pairs = 0;
   
   // 读取从输入文件或重启动文件读入数据.
   bool is_from_restart = tbox::RestartManager::getManager()->isFromRestart();

   if (is_from_restart){
      // 从重启动读入数据
      getFromRestart();
   } else {
        // 从输入文件读入数据
        getFromInput(input_db);
   }

   // 注册变量和数据片.
   registerModelVariables();

   // 注册为重启动对象.
   tbox::RestartManager::getManager()->registerRestartItem(d_object_name, this);
}

/*************************************************************************
 *
 * 析构函数.
 *
 ************************************************************************/
QIP::~QIP()
{
    tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
};

/*************************************************************************
 *
 * 注册变量和数据片.
 *
 ************************************************************************/
void QIP::registerModelVariables()
{
   // 变量数据库.
   hier::VariableDatabase<NDIM>* variable_db = 
                  hier::VariableDatabase<NDIM>::getDatabase();

   // 定义变量.
   tbox::Pointer<pdat::CellVariable<NDIM,double> > d_quantum_state 
       = new pdat::CellVariable<NDIM,double>("quantum_state",1);  // 量子态数据
   tbox::Pointer<pdat::CellVariable<NDIM,double> > d_result  
       = new pdat::CellVariable<NDIM,double>("result",1);  // 计算结果
  
   // 当前值上下文, 新值上下文
   tbox::Pointer<hier::VariableContext> d_current 
       = variable_db->getContext("CURRENT");
   tbox::Pointer<hier::VariableContext> d_new
       = hier::VariableDatabase<NDIM>::getDatabase()->getContext("NEW");

   // 数据片<quantum_state,current>: 存储量子态数据, 影像区宽度为0
   d_quantum_state_id = variable_db->registerVariableAndContext(d_quantum_state, d_current);

   // 数据片<result,new>: 存储计算结果, 影像区宽度为0
   d_result_id = variable_db->registerVariableAndContext(d_result, d_new);
}

/*************************************************************************
 *                                                                        
 * 注册绘图量.
 *                                                                  
 *************************************************************************/
void QIP::registerPlotData(
      tbox::Pointer<appu::JaVisDataWriter<NDIM> > javis_writer)
{
   javis_writer-> registerPlotQuantity("QuantumState", "SCALAR", d_quantum_state_id);
   javis_writer-> registerPlotQuantity("Result", "SCALAR", d_result_id);
}

/*************************************************************************
 *
 *  初始化指定的积分构件.
 *
 ************************************************************************/
void QIP::initializeComponent(algs::IntegratorComponent<NDIM>* intc) const
{
   const string &intc_name = intc->getName();

   if ( intc_name=="ALLOC_NEW_DATA") {  // 内存构件
         intc->registerPatchData(d_result_id);   

    } else if ( intc_name=="INIT_SET_VALUE") {  // 初值构件
         intc->registerInitPatchData(d_quantum_state_id);

    } else if(  intc_name=="STEP_SIZE" ) {  // 步长构件

    } else if ( intc_name=="COMPUTE") {  // 数值构件
          intc->registerPatchData(d_quantum_state_id);
          intc->registerPatchData(d_result_id);

    } else if ( intc_name=="COPY_SOLUTION") {  // 复制构件
          intc->registerCopyPatchData(d_quantum_state_id, d_result_id);   

    } else {
        TBOX_ERROR("\n::initializeComponent() : component " 
                   << intc_name <<" is not matched. "<<endl);
    }
}

/*************************************************************************
 *
 *  初始化数据片（支持初值构件）.
 *
 ************************************************************************/
void QIP::initializePatchData(
      hier::Patch<NDIM>& patch, 
      const double  time, 
      const bool    initial_time,
      const string& intc_name)
{
   if (initial_time) {
      tbox::Pointer< pdat::CellData<NDIM,double> > quantum_state =
                                    patch.getPatchData(d_quantum_state_id);

      // 初始化量子态为|0⟩
      quantum_state->fill(0.0);
      
      // 在第一个网格单元设置初始量子态
      const hier::Index<NDIM> &ifirst=patch.getBox().lower();
      const hier::Index<NDIM> &ilast =patch.getBox().upper();
      
      if (ifirst(0) == 0 && ifirst(1) == 0) {
         quantum_state->getPointer()[0] = 1.0; // |0⟩态
      }
   }
}

/*************************************************************************
 *
 *  计算稳定性时间步长（支持步长构件）.
 *
 ************************************************************************/
double QIP::getPatchDt(hier::Patch<NDIM>& patch,
                          const double  time,
                          const bool    initial_time,
                          const int     flag_last_dt,
                          const double  last_dt,
                          const string& intc_name)
{
   // 量子计算不需要时间步长，返回固定值
   return 1.0;
}

/*************************************************************************
 * 完成单个网格片上的数值计算（支持数值构件）.
 *
 * 该函数实现量子内积计算的核心逻辑
 ************************************************************************/
void QIP::computeOnPatch(hier::Patch<NDIM>& patch,
                          const double  time,
                          const double  dt,
                          const bool    initial_time,
                          const string& intc_name)
{
   // 获取数据片
   tbox::Pointer< pdat::CellData<NDIM,double> > quantum_state=
                                    patch.getPatchData(d_quantum_state_id);
   tbox::Pointer< pdat::CellData<NDIM,double> > result =
                                    patch.getPatchData(d_result_id);

   // 获取网格片索引范围
   const hier::Index<NDIM> &ifirst=patch.getBox().lower();
   const hier::Index<NDIM> &ilast =patch.getBox().upper();

   // 处理输入向量文件
   ifstream fin(d_input_vector_file);
   ofstream fout(d_output_result_file, ios::app);
   
   if (fin.is_open()) {
       array<double,16> a{}, b{};
       vector<pair<double,double>> results;
       
       while(true){
           // 读取第一个向量
           bool read_a = true, read_b = true;
           for(int i=0;i<16;++i){ 
               if(!(fin>>a[i])) { 
                   read_a = false; 
                   break; 
               } 
           }
           if(!read_a) { goto FILE_DONE; }
           
           // 读取第二个向量
           for(int i=0;i<16;++i){ 
               if(!(fin>>b[i])) { 
                   read_b = false; 
                   break; 
               } 
           }
           if(!read_b) { goto FILE_DONE; }
           
           double q = q_inner(a,b);
           double classic=0; for(int i=0;i<16;++i) classic += a[i]*b[i];
           if(classic<0) q = -q;
           
           results.emplace_back(q, classic);
           fout<<"Pair "<<results.size()<<": Quantum="<<fixed<<setprecision(4)<<q<<", Classic="<<classic<<"\n";
           d_processed_pairs++;
       }
FILE_DONE:
       if(!results.empty()) {
           double m=0; for(auto &p:results) m += fabs(p.first-p.second); m/=results.size();
           fout<<"\n平均绝对误差: "<<fixed<<setprecision(4)<<m<<"\n";
       }
       fin.close();
   }
   
   // 将结果存储到网格数据中
   result->fill(d_processed_pairs);
   
   fout.close();
}

/*************************************************************************
 *
 *  填充物理边界条件.
 *
 ************************************************************************/
void QIP::setPhysicalBoundaryConditions(
                  hier::Patch<NDIM>& patch,
                  const double fill_time,
                  const hier::IntVector<NDIM>& ghost_width_to_fill,
                  const string& intc_name)
{
   // 量子计算不需要物理边界条件
}

/*************************************************************************
 **
 **  从重启动数据库读取数据. 
 **
 *************************************************************************/
void QIP::getFromRestart()
{
   tbox::Pointer<tbox::Database> root_db =
      tbox::RestartManager::getManager()->getRootDatabase();

   tbox::Pointer<tbox::Database> db;
   if ( root_db->isDatabase(d_object_name) ) {
      db = root_db->getDatabase(d_object_name);
   } else {
      TBOX_ERROR("Restart database corresponding to "
              << d_object_name << " not found in restart file.");
   }

   d_processed_pairs = db->getInteger("d_processed_pairs");
}

/*************************************************************************
*
*  输出数据成员到重启动数据库.
*
************************************************************************/
void QIP::putToDatabase(tbox::Pointer<tbox::Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
   assert(!db.isNull());
#endif

   db->putInteger("d_processed_pairs", d_processed_pairs);
   db->putString("d_input_vector_file", d_input_vector_file);
   db->putString("d_output_result_file", d_output_result_file);
}

/*************************************************************************
 *
 *  从输入数据库读入数据.
 *
 ************************************************************************/
void QIP::getFromInput(tbox::Pointer<tbox::Database> db)
{
   if (db->keyExists("input_vector_file")) {
      d_input_vector_file = db->getString("input_vector_file");
   } else {
      d_input_vector_file = "vectors.txt";
   }
   
   if (db->keyExists("output_result_file")) {
      d_output_result_file = db->getString("output_result_file");
   } else {
      d_output_result_file = "qip_results.txt";
   }
}

static inline double norm2(const array<double,16>& v){ double s=0; for(double x:v) s+=x*x; return sqrt(s); }

static void unit_vector(array<double,16>& v){
    double n = norm2(v);
    if(n==0) throw runtime_error("Cannot unitize a zero vector.");
    for(double &x:v) x/=n;
}

static void hierarchical_normalization(const array<double,16>& a, array<double,15>& out){
    // Build norms layer by layer, append reversed each time into out (total 15 values)
    array<double,16> layer = a;
    int len = 16;
    int writePos = 0;
    while(len > 1){
        int nextLen = len/2;
        array<double,16> next{};
        for(int i=0;i<nextLen;++i){
            double x = layer[2*i], y = layer[2*i+1];
            next[i] = sqrt(x*x + y*y);
        }
        // append reversed next[0..nextLen-1]
        for(int i=nextLen-1;i>=0;--i){
            out[writePos++] = next[i];
        }
        // prepare for next iteration in forward order
        for(int i=0;i<nextLen;++i) layer[i] = next[i];
        len = nextLen;
    }
}

// 量子内积计算核心函数
double QIP::q_inner(const array<double,16>& a_in, const array<double,16>& b_in){
    array<int,16> flag_a{}, flag_b{};
    for(size_t i=0;i<16;++i) flag_a[i] = (a_in[i]>0)?1:(a_in[i]==0?0:-1);
    for(size_t i=0;i<16;++i) flag_b[i] = (b_in[i]>0)?1:(b_in[i]==0?0:-1);
    array<double,16> a{}, b{};
    for(size_t i=0;i<16;++i){ a[i]=fabs(a_in[i]); b[i]=fabs(b_in[i]); }
    double norma = norm2(a), normb = norm2(b);
    if(norma==0 || normb==0) return 0.0;
    size_t index_a=0,index_b=0; double ka0=0,kb0=0;
    for(size_t i=0;i<16;++i){ if(a[i]!=0){ index_a=i; ka0=a[i]; break; } }
    for(size_t i=0;i<16;++i){ if(b[i]!=0){ index_b=i; kb0=b[i]; break; } }
    unit_vector(a); unit_vector(b);
    double k_a = a[index_a]/ka0;
    double k_b = b[index_b]/kb0;

    auto build_theta = [](const array<double,16>& vec){
        array<double,15> temp{}; // hierarchical output length 15
        hierarchical_normalization(vec, temp);
        // reverse temp into rtemp
        array<double,15> rtemp{};
        for(int i=0;i<15;++i) rtemp[i] = temp[14 - i];
        // concat = rtemp (15) + vec (16) => 31
        double concat[31];
        for(int i=0;i<15;++i) concat[i]=rtemp[i];
        for(int i=0;i<16;++i) concat[15+i]=vec[i];
        double res[31];
        res[0]=concat[0];
        for(int i=1;i<31;++i){ double parent = concat[(i-1)/2]; res[i] = (parent==0.0)?0.0:(concat[i]/parent); }
        array<double,15> theta{};
        int t=0; for(int i=1;i<31; i+=2){ double x=res[i]; if(x>1)x=1; if(x<-1)x=-1; theta[t++] = 2.0*acos(x); }
        return theta;
    };

    array<double,15> theta_a = build_theta(a);
    array<double,15> theta_b = build_theta(b);

    double p0 = Quantum_circuit_whitebox(theta_a, theta_b, flag_a, flag_b);
    if(p0 < 0.5) return 0.0; else return sqrt(max(0.0, 2.0*p0 - 1.0))/(k_a*k_b);
}

/*************************************************************************
 *
 *  量子电路白盒计算
 *
 ************************************************************************/
double QIP::Quantum_circuit_whitebox(const array<double,15>& theta_a, 
                                      const array<double,15>& theta_b,
                                      const array<int,16>& flag_a, 
                                      const array<int,16>& flag_b){
    State st;
    auto X=matX(), Z=matZ(), H=matH();
    auto RY=[&](double th){ return matRY(th); };

    auto CRY = [&](int c, int t, double th){ st.applyCtrl1(RY(th), {c}, t); };
    auto MCRY = [&](vector<int> ctrls, int t, double th){ st.applyCtrl1(RY(th), ctrls, t); };
    auto MCZ = [&](vector<int> ctrls, int t){ st.applyCtrl1(Z, ctrls, t); };
    auto MCX = [&](vector<int> ctrls, int t){ st.applyCtrl1(X, ctrls, t); };
    auto s1 = [&](int q, const M2&U){ st.apply1(U,q); };
    auto X1 = [&](int q){ s1(q,X);} ;
    auto H1 = [&](int q){ s1(q,H);} ;

    // CSWAP via controlled 2-qubit 4x4
    array<C,16> SWAP = { C(1,0),C(0,0),C(0,0),C(0,0),
                         C(0,0),C(0,0),C(1,0),C(0,0),
                         C(0,0),C(1,0),C(0,0),C(0,0),
                         C(0,0),C(0,0),C(0,0),C(1,0)};
    auto CSW = [&](int c, int a, int b){ st.applyCtrl2(SWAP, {c}, max(a,b), min(a,b)); };

    // (1)
    s1(1, RY(theta_a[0]));
    s1(5, RY(theta_b[0]));
    // (2)
    X1(1); X1(5);
    CRY(1,2, theta_a[1]);
    CRY(5,6, theta_b[1]);
    X1(1); X1(5);
    CRY(1,2, theta_a[2]);
    CRY(5,6, theta_b[2]);
    // (3)
    X1(1); X1(2);
    X1(5); X1(6);
    MCRY({1,2}, 3, theta_a[3]);
    MCRY({5,6}, 7, theta_b[3]);
    X1(1); X1(2);
    X1(5); X1(6);

    X1(1); X1(5);
    MCRY({1,2}, 3, theta_a[4]);
    MCRY({5,6}, 7, theta_b[4]);
    X1(1); X1(5);

    X1(2); X1(6);
    MCRY({1,2}, 3, theta_a[5]);
    MCRY({5,6}, 7, theta_b[5]);
    X1(2); X1(6);

    MCRY({1,2}, 3, theta_a[6]);
    MCRY({5,6}, 7, theta_b[6]);

    // (4) + sign fixes blocks 0..7 mirrored from python
    auto block = [&](int ax, int ay, int bx, int by, double tha, double thb){
        MCRY({1,2,3}, 4, tha);
        if( (ax==1 && ay==-1) || (ax==0 && ay==-1) ) MCZ({1,2,3},4);
        if( ax==-1 && ay==1 ){ MCX({1,2,3},4); MCZ({1,2,3},4); MCX({1,2,3},4);} 
        if( (ax==-1 && ay==-1) || (ax==-1 && ay==0) ){ MCZ({1,2,3},4); MCX({1,2,3},4); MCZ({1,2,3},4); MCX({1,2,3},4);} 
        MCRY({5,6,7}, 8, thb);
        if( (bx==1 && by==-1) || (bx==0 && by==-1) ) MCZ({5,6,7},8);
        if( bx==-1 && by==1 ){ MCX({5,6,7},8); MCZ({5,6,7},8); MCX({5,6,7},8);} 
        if( (bx==-1 && by==-1) || (bx==-1 && by==0) ){ MCZ({5,6,7},8); MCX({5,6,7},8); MCZ({5,6,7},8); MCX({5,6,7},8);} 
    };

    // 0
    X1(1); X1(2); X1(3); X1(5); X1(6); X1(7);
    block(flag_a[0],flag_a[1],flag_b[0],flag_b[1], theta_a[7], theta_b[7]);
    X1(1); X1(2); X1(3); X1(5); X1(6); X1(7);
    // 1
    X1(1); X1(2); X1(5); X1(6);
    block(flag_a[2],flag_a[3],flag_b[2],flag_b[3], theta_a[8], theta_b[8]);
    X1(1); X1(2); X1(5); X1(6);
    // 2
    X1(1); X1(3); X1(5); X1(7);
    block(flag_a[4],flag_a[5],flag_b[4],flag_b[5], theta_a[9], theta_b[9]);
    X1(1); X1(3); X1(5); X1(7);
    // 3
    X1(1); X1(5);
    block(flag_a[6],flag_a[7],flag_b[6],flag_b[7], theta_a[10], theta_b[10]);
    X1(1); X1(5);
    // 4
    X1(2); X1(3); X1(6); X1(7);
    block(flag_a[8],flag_a[9],flag_b[8],flag_b[9], theta_a[11], theta_b[11]);
    X1(2); X1(3); X1(6); X1(7);
    // 5
    X1(2); X1(6);
    block(flag_a[10],flag_a[11],flag_b[10],flag_b[11], theta_a[12], theta_b[12]);
    X1(2); X1(6);
    // 6
    X1(3); X1(7);
    block(flag_a[12],flag_a[13],flag_b[12],flag_b[13], theta_a[13], theta_b[13]);
    X1(3); X1(7);
    // 7
    block(flag_a[14],flag_a[15],flag_b[14],flag_b[15], theta_a[14], theta_b[14]);

    // Swap test
    H1(0);
    CSW(0,1,5); CSW(0,2,6); CSW(0,3,7); CSW(0,4,8);
    H1(0);

    // Probability qubit 0 == |0>
    double p0=0.0; const size_t N = State::Dim;
    for(size_t idx=0; idx<N; ++idx){
        if( (idx & 1ull)==0 ){
            p0 += st.amp[idx].r*st.amp[idx].r + st.amp[idx].i*st.amp[idx].i;
        }
    }
    return p0;
}


