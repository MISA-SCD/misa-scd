<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<script type="text/x-mathjax-config">
    MathJax.Hub.Config({ tex2jax: {inlineMath: [['$', '$']]}, messageStyle: "none" });
</script>
# MISA-SCD v1.0用户文档

## 1. 概述
### 1.1. MISA-SCD简介
MISA-SCD v1.0是一款并行空间分辨随机团簇动力学（Spatially Resolved Stochastic Cluster Dynamics，SRSCD）模拟软件，主要用于模拟纯Fe、FeCu合金材料在辐照、热老化及退火时的缺陷团簇长时间演化行为，获得缺陷团簇的成分、尺寸、数密度分布等信息，用以预测材料的宏观性能变化，后续将增加对其他材料（如钨机材料）的模拟支持。其中辐照条件包括电子辐照、中子辐照。缺陷团簇的类型包括空位（V）、自间隙原子（SIA）、溶质原子（Cu）、空位团簇、自间隙团簇、溶质团簇。

### 1.2. MISA-SCD文件结构

#### `inputs`：
* `cascades`：级联缺陷文件，包含用于中子辐照模拟所需的级联缺陷的个数及空间分布信息。
* `defectsAttributes`：缺陷特征能量参数文件，包含用于各种模拟条件所需的缺陷特征能量参数。
* `mesh`：网格文件：包含用于模拟网格配置信息（网格为立方体），包括box的边界条件，网格尺寸，x、y、z方向的网格数量。
#### `src`：
  包含所有源代码，其中src_MISASCD.f90是主程序。
#### `tests`：
  包含各算例，其中,每个算例配备一个parameters.txt，用于模拟的输入参数配置。  
  RPV_FeCu_Cascade：FeCu合金，中子辐照模拟  
  RPV_FeCu_FrenkelPairs：FeCu合金，电子辐照模拟  
  SANS_Fe_Cascade：纯Fe，中子辐照模拟  
  SANS_Fe0.1Cu_Cascade：Fe-0.1%at.Cu合金，中子辐照模拟  
  SANS_Fe0.3Cu_Cascade：Fe-0.3%at.Cu合金，中子辐照模拟  
  SRSCD_Fe_Cascade：纯Fe，中子辐照模拟  
  
## 2. 编译

### 2.1. 准备环境
* MPI环境
* 支持Fortran90的Fortran编译器

### 2.2. 编译
使用makefile编译， makefile位于src目录下。

```bash
$ module load MPI_MODULE  #加载mpi环境
$ cd $MISA_SCD_PATH/src/   #$MISA_SCD_PATH为misa-scd的路径
$ make
```
编译完成后，在src目录下会生成二进制可执行文件`misascd`。

## 3. 运行
### 3.1. 准备输入配置文件
MISA-SCD的输入配置文件包括三类：**基准参数配置文件**、**网格参数配置文件**、**缺陷特征能量参数配置文件**。其中，对于中子辐照模拟，还需包括**级联缺陷文件**，各配置文件说明如下
* **基准参数配置文件** 位于各算例目录下，用户可根据模拟需求进行修改。
* **网格参数配置文件** 位于inputs/mesh/目录下，用户可根据模拟需要的体系大小进行修改。在使用时，需在基准参数配置文件中，对应参数下方给出网格参数配置文件的相对路径。
* **缺陷特征能量参数文件** 位于inputs/defectsAttributes/目录下，该目录下已提供若干缺陷特征能量参数文件，基本能够满足纯Fe、FeCu体系在辐照、热老化及退火条件下的模拟需求。在使用时，只需在基准参数配置文件中，对应参数下方给出缺陷特征能量参数文件的相对路径即可。如需修改能量参数请进入对应的文件中修改，关于文件的说明详见文件头部信息。用户也可以创建自己的参数文件，格式参考inputs/defectsAttributes/目录中给出的例子。
* **级联缺陷文件** 位于inputs/cascades/目录下，该文件中的信息一般由分子动力学模拟级联碰撞过程获得。在使用时，需在基准参数配置文件中，对应参数下方给出级联缺陷文件的相对路径。在inputs/cascades/目录下已提供一个20KeV PKA产生的级联缺陷文件。用户也可以创建自己的级联缺陷文件，格式参考提供的例子。
  
关于配置文件中各参数的相关说明请参考**输入配置项说明**章节。
### 3.2. 运行MISA-SCD
> 在相应算例的文件夹下运行可在执行文件misascd.  
> 注意：进程数不能小于网格数，要保证每个进程至少一个网格。运行不同规模的模拟之前，先检查inputs/meshes中对应的网格配置文件中的网格数是否满足此条件。

MISA-SCD可通过`mpirun`命令运行，下面的示例中，使用8个MPI进程，parameters.txt为基准参数配置文件。
```bash
$ cd $MISA_SCD_PATH/tests/example_test_filename/
$ mpirun -n 8 ../../src/misascd parameters.txt   #实际使用时../../src/misascd替换为可执行文件misa的实际路径
```
也可通过`sbatch`提交脚本来运行MISA-SCD，要求运行环境已安装配置`SLURM`。

```bash
$ cd $MISA_SCD_PATH/tests/example_test_filename/
$ sbatch misascd.sh
```

脚本misascd.sh如下：
```bash
$ cat misascd.sh
$ #!bin/sh
$ #SBATCH -J jobName
$ #SBATCH -O %j.out
$ #SBATCH -p test
$ #SBATCH -t 10:00:00
$ #SBATCH -N 1
$ #SBATCH -n 32

mpirun -n 32 ../../src/misascd parameters.txt   #实际使用时../../src/misascd替换为可执行文件misa的实际路径
```

## 4. 用MISA-SCD模拟缺陷演化过程
### 4.1. 基本原理
MISA-SCD基于空间分辨随机团簇动力学（Spatially Resolved Stochastic Cluster Dynamics，SRSCD）方法建模，用于模拟点缺陷及其团簇在辐照、热老化、退火条件下的产生、扩散、聚集、分解、湮灭等动力学行为，获取不同时刻下缺陷的数密度和尺寸信息。在MISA-SCD中，整个模拟区域被划分为均匀的小体积元（即网格，如图 4-1 所示），体积元内的缺陷均匀分布，可动缺陷可在体积元之间扩散。缺陷之间的相互作用被视为一个个反应（类似化学反应），主要包括Insertion（注入Frenkel对或级联缺陷）、Dissocation(缺陷团簇分解出一个点缺陷而变小)、SinkRemoval（被材料内固有缺陷捕获而消失）、Clustering（缺陷团簇吸收点缺陷/小团簇而变大）、Diffusion（可动缺从一个体积元扩散到另一个体积元），共五类反应，如图 4-2 所示。其中，Insertion反应仅发生在辐照条件下。

在并行方面，MISA-SCD采用同步并行KMC算法（Synchronous Parallel Kinetic Monte Carlo，SPKMC）进行进程间的时间同步。在每个时间步内，利用KMC直接法在当前进程内选择反应并发生，具体步骤如下：
1. 生成一个$(0,1]$的随机数，为所有进程选择一个时间增量$\delta t$；
2. 生成一个随机数，在当前进程内随机选择一个体积元；
3. 生成一个$(0,1]$的随机数，在所选体积元内随机选择一个反应；
4. 根据选择的反映更新相应的缺陷，及缺陷涉及的反应；
5. 体系时间前进$\delta t$。

<div  align="center">    
<img src="./MISA-SCDBox.png" width = "60%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>
<center>图 4-1 MISA-SCD空间区域划分示意图</center>

<div  align="center">    
<img src="./reactions.jpg" width = "50%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>
<center>图 4-2 MISA-SCD中各类反应示意图</center>

### 4.2. 术语解释
* 点缺陷：只含有一个空位（自间隙原子）的缺陷。
* 团簇：含有多个点缺陷的缺陷。
* 缺陷尺寸：含有点缺陷的个数即为该缺陷的尺寸。
* Frenkel对：包括一个空位和一个自间隙原子，电子辐照下，当发生Insertion反应时，向体积元内放入一个Frenkel对。
* 级联缺陷：由级联碰撞产生，包括多个缺陷，有点缺陷，还有团簇。
* KMC：动力学蒙特卡洛
* DPA：Displavement Per Atom，用于衡量辐照剂量。

下面以**FeCu体系中子辐照模拟**为例，说明MISA-SCD的使用过程。

### 4.3. 输入配置
#### 配置缺陷特征能量
FeCu体系中子辐照模拟的缺陷特征能量参数文件为FeCu_Defects_Cas.txt，该文件中主要包含：缺陷中所含点缺陷的最大类型数、点缺陷的类型及其形成能（$E_f$）、点缺陷及小团簇的类型及其迁移能（$E_m$）和扩散前置因子（$D_0$）（用于计算该缺陷的扩散率（$D$））、用于计算大团簇扩散率的信息、小团簇的类型及其结合能（$E_b$）、用于计算大团簇结合能的信息、允许的反应类型及参与反应的缺陷类型信息。下面逐一说明：
1. **缺陷类型：**
<div  align="center">    
<img src="./species.png" width = "30%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>
默认为4，即MISA-SCD中用一个长度为4的数组存储缺陷类型，数组第一个元素表示该缺陷中所含Cu原子个数，第二个元素表该缺陷中所含空位个数，第三个元素中表示该团簇中所含自间隙原子的个数（该值大于0表示该缺陷为可动的自间隙团簇），第四个元素同样表示该团徐中所含自间隙原子的个数（该值大于0表示该缺陷为不可动的自间隙团簇）。

2. **点缺陷的类型及其形成能（$E_f$）：**
<div  align="center">    
<img src="./Ef.png" width = "60%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

* 57行：`formationEnergies`为标识，表示下面开始是缺陷的形成能信息。  
* 59~60行：`numSingle`的值为3，表示有3个形成能。  
* 62~63：缺陷类型及其对应的形成能，`1 0 0 0`表示该缺陷为Cu原子，其对应的形成能为1.77eV。

3. **点缺陷及小团簇的类型及其迁移能（$E_m$）和扩散前置因子（$D_0$）：**
<div  align="center">    
<img src="./Em.png" width = "40%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

* 69行：`diffusionPrefactors`为标识，表示下面开始是缺陷的迁移能和扩散前置因子信息。  
* 71~72行：`numSingle`的值为3，表示有3个迁移能和扩散前置因子。  
* 74~75：缺陷类型及其对应的扩散前置因子和迁移能，`1 0 0 0`表示该缺陷为Cu原子，其对应的扩散前置因子为$6.3\times10^{13} {nm}^2/s$，迁移能为2.29eV。

4. **用于计算大团簇扩散率的信息：**
<div  align="center">    
<img src="./EmFunc.png" width = "70%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>
在MISA-SCD中，大团簇的扩散率由经验公式计算得出，不同的缺陷，计算公式不同。  

* 81~82行：`numFunction`的值为4，表示有4个。  
* 84行：缺陷类型的基本表示，表示该类缺陷是仅含有Cu原子的团簇。  
* 85~86行：该类缺陷的最小尺寸和最大尺寸，-1表示无穷大。  
* 87行：计算扩散率的公式标识，对于不同的缺陷类型，有不同的计算公式。  
* 88行：计算该类缺陷的扩散率所需的参数个数，这里为0个，表示该类缺陷的扩散率为0 ，即尺寸大于等于2的Cu团簇不可动。若参数个数大于0个，则在下一行列出参数，并以空格隔开。

5. **小团簇的类型及其结合能（$E_b$）：**
<div  align="center">    
<img src="./Eb.png" width = "50%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

结合能是用来计算分解反应的反应速率的。
* 108行：`diffusionPrefactors`为标识，表示下面开始是缺陷的迁移能和扩散前置因子信息。  
* 110~111行：`numSingle`的值为46，表示46个结合能。  
* 113~114：缺陷类型及其对应的结合能，`2 0 0 0`表示该缺陷为含2个Cu原子Cu团簇，其分解出去1个Cu原子（即`1 0 0 0`）所需的结合能是0.19eV。

6. **用于计算大团簇结合能的信息：**
<div  align="center">    
<img src="./EbFunc.png" width = "50%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>
在MISA-SCD中，大团簇的结合能由经验公式计算得出，不同的缺陷，计算公式不同。  

* 206~207行：`numFunction`的值为11，表示有11个。  
* 209行：缺陷类型的基本表示，第一个`1 0 0 0`表示该类缺陷是仅含有Cu原子的团簇，第二个`1 0 0 0`表示该Cu团簇分解出去一个Cu原子。 
* 210~211行：该类缺陷的最小尺寸和最大尺寸，-1表示无穷大。  
* 212行：计算结合能的公式标识，对于不同的缺陷类型，有不同的计算公式。  
* 213行：计算该类缺陷的结合能所需的参数个数，这里为3个。
* 214行：依次列出所需的3个参数的值。

7. **允许的反应类型及参与反应的缺陷类型信息**
   
MISA-SCD中，不同反应的反应速率由不同的公式计算得出，具体采用哪个公式，由该文件中的`fType`给出。
* Dissociation反应：
参与Dissociation反应的反应物为1个，产物为2个（一个为点缺陷，另一个为产物为反应物分解出1个点缺陷后转变成为的缺陷）。
<div  align="center">    
<img src="./Dissoc.png" width = "50%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

289行：反应标识。  
290行：Dissociation反应的个数。   
292行：反应物及产物的基本类型，第一个`1 0 0 0`表示参与Dissociation反应的反应物的基本类型（这里为Cu团簇）；第二个`1 0 0 0`表示参与Dissociation反应的其中一个产物的类型（这里为Cu原子）。这里表示Cu团簇发生Dissociation 反应，分解出去一个Cu原子，则该Cu团簇尺寸减1（即第二个产物）。  
293~294行：反应物的最小尺寸和最大尺寸，-1表示无穷大。  
295行：计算反应速率的公式标识。

* Diffusion反应：
参与Diffusion反应的反应物和产物都为1个，缺陷类型不变。
<div  align="center">    
<img src="./diff.png" width = "50%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

317行：反应标识。  
318行：Diffusion反应的个数。  
320行：反应物的基本类型，第一个`1 0 0 0`表示反应物的基本类型；第二个`1 0 0 0`表示产物的基本类型。  
321~322行：反应物的最小尺寸和最大尺寸。  
323行：计算反应速率的公式标识。

* SinkRemoval反应：
参与SinkRemoval反应的反应物为1个，产物为0个。
<div  align="center">    
<img src="./sink.png" width = "50%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

335行：反应标识。  
336行：SinkRemoval反应的个数。  
338行：反应物的基本类型，这里标识空位团簇。 
339~340行：反应物的最小尺寸和最大尺寸。  
341行：计算反应速率的公式标识。

* Clustering反应：
参与Clustering反应的反应物为2个，产物为0个/1个/2个，具体视反应物的类型而定。
<div  align="center">    
<img src="./clu.png" width = "50%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

356行：反应标识。  
357行：SinkRemoval反应的个数。  
359行：两个反应物的基本类型。  
360~361行：两个反应物的最小尺寸和最大尺寸，-1标识无穷大。   
362行：计算反应速率的公式标识。

### 配置级联缺陷文件
该参数文件仅中子辐照模拟时用到。一个文件中可以包含多组级联缺陷信息，一般是同样PKA能量下的级联碰撞模拟获得的。
<div  align="center">    
<img src="./cas1.png" width = "100%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

1行：平均移位原子数量。  
3行：该文件中由9组级联缺陷。  
下面依次列出各组级联缺陷：
<div  align="center">    
<img src="./cas2.png" width = "60%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

5行：该组级联缺陷所包含的缺陷个数。  
6行：该组级联缺陷所包含的移位原子个数。 
下面依次列出该组中的级联缺陷： 
7~8行：缺陷类型，及其偏移坐标（x，y，z方向），即相对于级联中心的偏移坐标。当一组级联缺陷被注入某网格时，其级联中心即为该网格的中心。

#### 配置基准输入
该文件中，每个参数包含两行数据，第一行为该参数的名称，第二行为该参数的值，两个参数之间以空行隔开。主要配置以下参数：
1. 控制参数：包括缺陷特征能量参数文件、网格文件、级联缺陷文件的相对路径，缺陷注入格式、是否考虑晶界、是否仅点缺陷可动等；
2. 模拟参数：包括模拟温度，合金含量、dpa速率、总dpa、位错密度、晶格常数、原子体积、模拟次数等；
3. 输出控制参数：包括是否输出totdat.txt文件、参与输出统计的团簇最小尺寸；
4. 细网格参数：包括细网格的边长，x、y、z方向的细网格的数量。该项仅在中子辐照、且，缺陷注入格式为`adaptive`时需配置。
该文件中各配置项的详细说明参见**输入配置项说明**章节。

#### 配置网格信息
该文件中各配置项的详细说明参见**输入配置项说明**章节。

### 4.4. 输出结果分析
输出文件包括两个：屏幕输出和totdat X.out，其中“X”表示重复模拟时，第X次模拟的输出文件。对于一次模拟过程，中间时刻的结果和最终时刻的结果都写在同一个totdat X.out文件中，两次输出时刻的结果之间，用空行隔开，如下所示。

<div  align="center">    
<img src="./output.png" width = "90%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>

包括三部分信息：  
2~5行：输出时刻信息；  
6~10行：缺陷信息；  
11~22行：统计信息。

其中缺陷信息的前4列为缺陷类型，最后一列为该类缺陷的数量（整体模拟空间中的）。统计信息中主要统计四类缺陷：S（纯Cu图团簇、含有Cu和空位的团簇）、Void（空位团簇）、Loop（自间隙团簇）、SV（含有Cu和空位的团簇）。

对于MISA-SCD的结果进行分析，一般是获取某类缺陷团簇的总数密度、平均半径随时间（DPA）变化的情况。可采用两种方式进行统计分析，一种是在基准输入参数配置文件中设置**输出参数控制项**，由MISA-SCD的后处理程序自动统计不同输出时刻下某类团簇的总数密度及平均尺寸，并写入totdat X.out文件中（写在缺陷信息下方）；另一种是根据totdat X.out中不同输出时刻输出的缺陷团簇信息（缺陷类型及其数量），统计某类缺陷总数密度及平均尺寸。


获取不同输出时刻下某类缺陷团簇的总数密度及平均尺寸数据后，可采用orgin软件进行绘图展示，下图给出了一个示例。

<div  align="center">    
<img src="./CuNumberDensity-Cas.jpg" width = "50%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>
<center>图 4-3 Cu团簇数密度随时间和DPA的变化</center>

<div  align="center">    
<img src="./CuAverageRadius-Cas.jpg" width = "50%" alt="MISA-SCD模拟区域分解示意图" align=center />
</div>
<center>图 4-4 Cu团簇平均半径随时间和DPA的变化</center>

## 5. 输入配置项说明
这里仅给出**基准参数配置文件**及**网格参数配置文件**中各配置项的说明，对于**缺陷特征能量参数配置文件**
及**级联缺陷文件**的配置说明请参见章节**4. 示例**。
### 5.1. 基准参数配置文件
例如：parameters.txt，其中的信息分为三部分：**控制参数**、**模拟参数**、**输出控制参数**、**细网格参数**，各部分参数项说明如下：

**-----------控制参数项-----------**

***defectFile***  
类型：字符串  
说明：必选项。缺陷特征能量参数文件的相对路径，即inputs/defectsAttributes/中的对应文件。

***meshFile***  
类型：字符串  
说明：必选项。网格文件的相对路径，即inputs/meshes/中的对应文件。

***implantScheme***  
类型：字符串  
说明：可选项。若模拟辐照条件，则该字段必须设置；若模拟热老化及退火，则无须设置该字段。该字段表示初始缺陷置入box中的方式，包括`MonteCarlo` 或 `explicit` 两种方式，默认为 `MonteCarlo`。`MonteCarlo` 表示随机置入，置入的时刻及置入哪个网格，都是随机的；`explicit` 表示每隔指定的模拟时长时置入，置入哪个网格则是随机的。注意：`MonteCarlo` 可用于各种辐照条件的模拟，`explicit` 仅用于中子辐照的模拟。对于电子辐照模拟，置入的初始缺陷为Frenkel对（一个空位和一个自间隙原子）；对于中子辐照模拟，置入的初始缺陷为级联缺陷（即**级联缺陷文件**中的信息，包括若干缺陷团簇，有空位、自间隙原子、空位团簇、自间隙团簇）。

***cascadeFile***  
类型：字符串  
说明：可选项。级联缺陷文件的相对路径。即inputs/cascades/中的对应文件。若为中子辐照模拟，则必须提供该文件。

***meshingType***  
类型：字符串  
说明：必选项。是否使用“自适应网格”，`adaptive` 或 `nonAdaptive` 两种方式，默认为 `nonAdaptive`。`adaptive` 表示使用“自适应网格”，仅用于中子辐照模拟；`nonAdaptive` 表示不使用“自适应网格”，用于各种辐照条件的模拟。若该字段设置为`adaptive`，则需设置下述**细网格参数项**中的各参数。

***grainBoundaries***  
类型：字符串  
说明：必选项。是否考虑晶界，`yes` 或 `no`。

***pointDefect***  
类型：字符串  
说明：必选项。是否仅点缺陷可动，`yes` 或 `no`。该字段与**缺陷特征能量参数**文件中的信息关联。详见**缺陷特征能量参数配置文件**章节中的说明。

**-----------模拟参数项-----------**

***temperature***  
类型：double precision  
说明：必选项。模拟温度，单位：K。一般不低于室温。

***CuContent***  
类型：double precision  
说明：必选项。FeCu合金中的Cu含量，单位：无。对于纯Fe体系，需设置为0；对于FeCu体系，需设置大于0的值。

***numVac***  
类型：integer  
说明：必选项。初始时体系中的空位个数，单位：无。默认为0

***dpaRate***  
类型：double precision  
说明：必选项。DPA速率，单位：dpa/s。对于辐照条件，需设置大于0的值；对于热老化及退火模拟，需设置为0。

***totalDPA***  
类型：double precision  
说明：必选项。总的辐照剂量，单位：dpa。对于辐照模拟，需设置大于0的值；对于热老化及退退火模拟，需设置为0.

***agingTime***
类型：double precision
说明：可选项。热老化模拟时长，单位：s（秒）。对于热老化模拟，需设置大于0的值；对于辐照及退火模拟，无须设置该值（任意值均可，或者直接删除该字段，但其值不能为空）。


***firr***  
类型：double precision  
说明：可选项。辐照增强因子，单位：无。对于FeCu体系模拟，需要设置大于等于1的值；对于纯Fe体系，无须设置该值（任意值均可，或者直接删除该字段，但其值不能为空）。

***atomSize***  
类型：double precision  
说明：必选项。Fe原子体积，单位：$nm^{-3}$。需设置大于0 的值。

***lattice***  
类型：double precision  
说明：必选项。晶格常数，单位：$nm$。需设置大于0 的值。对于Fe基材料，其值默认为0.01178$nm^3$。

***burgers***  
类型：double precision  
说明：必选项。伯格斯矢量，单位：$nm$。需设置大于0 的值。对于Fe基材料，其值默认为0.2867nm。

***reactionRadius***  
类型：double precision  
说明：必选项。复合半径，单位：$nm$。需设置大于0 的值。

***annealTemp***  
类型：double precision  
说明：可选项。退火温度，单位：K。对于退火模拟，需设置该字段的值，一般大于等于室温；对于辐照及热老化模拟，无须设置该值（任意值均可，或者直接删除该字段，但其值不能为空）。

***annealSteps***  
类型：integer  
说明：可选项。退火步数，单位：无。对于退火模拟，需设置大于等于1的值；对于辐照及热老化模拟，无须设置该值（任意值均可，或者直接删除该字段，但其值不能为空）。

***annealTime***  
类型：double precision  
说明：可选项。退火时间，单位：s（秒）。对于退火模拟，需设置大于0的值；对于辐照及热老化模拟，无须设置该值（任意值均可，或者直接删除该字段，但其值不能为空）。

***annealType***  
类型：字符串  
说明：可选项。退火时的温度变化类型，`add` 或 `mult` 两种方式。`add` 表示以等差数列的方式调整退火温度；`mult` 表示以等比数列的方式调整退火温度。调整的总次数等于 `annealSteps`的值。对于退火模拟，需设置该字段的值；对于辐照及热老化模拟，无须设置该值（`add` 或 `mult` 均可，或者直接删除该字段，但其值不能为空）。
若为`add`，则annealTemp = annealTemp + annealTempInc；若为`mult`，则annealTemp = annealTemp * annealTempInc。

***annealTempInc***  
类型：double precision  
说明：可选项。退火时的温度调整的增量，单位：K。对于退火模拟，需设置该字段的值（正、负、0均可）；对于辐照及热老化模拟，无须设置该值（任意值均可，或者直接删除该字段，但其值不能为空）。

***grainSize***  
类型：double precision  
说明：必选项。晶粒尺寸，单位：nm。需设置大于0的值。

***dislocDensity***  
类型：double precision  
说明：必选项。位错密度，单位：$nm^{-2}$。需设置大于0 的值。

***impurityConc***  
类型：double precision  
说明：可选项。杂质浓度，单位：${atom}^{-1}$。默认值为0。当前版本的程序暂不支持杂质加入的模拟。

***max3DInt***  
类型：integer  
说明：必选项。自间隙团簇的尺寸小于等于该值，则视为球状团簇；否则视为环状，单位：无。需设置大于等于1的值。该字段与**缺陷特征能量参数文件**中的信息关联，详见**缺陷特征能量参数配置**章节的相关说明。

***cascadeVolume***  
类型：double precision  
说明：可选项。级联范围，单位：$nm^{3}$。对于中子辐照模拟，要保证每个网格的体积大于等于该值；对于热老化及退火模拟，无须设置该值（任意值均可，或者直接删除该字段，但其值不能为空）。

***numSims***  
类型：integer  
说明：可选项。重复模拟的次数，单位：无。默认为1.

**-----------输出控制参数项-----------**

***totdatToggle***  
类型：字符串  
说明：可选项。是否输出缺陷信息统计文件totdat.out，`yes` 或 `no`。

***minSCluster***  
类型：integr  
说明：可选项。所含Cu原子个数大于该值的Cu团簇参与统计。若**totdatToggle**设置为`yes`，则必须设置该值。

***minVoid***  
类型：integr  
说明：可选项。所含空位个数大于该值的空位团簇参与统计。若**totdatToggle**设置为`yes`，则必须设置该值。

***minLoop***  
类型：integr  
说明：可选项。所含自间隙原子个数大于该值的自间隙团簇参与统计。若**totdatToggle**设置为`yes`，则必须设置该值。

***minSV***  
类型：integr  
说明：可选项。所含Cu原子个数及空位个数之和大于该值的Cu_空位复合团簇参与统计。若**totdatToggle**设置为`yea`，则必须设置该值。

**-----------细网格参数项-----------**

***fineLength***  
类型：double precision  
说明：可选项。细网格的尺寸，单位：nm。当 **meshingType** 的值为 `adaptive` 时，需要设置该值。

***numxFine***  
类型：integer  
说明：可选项。x方向的细网格个数，单位：无。当 **meshingType** 的值为 `adaptive` 时，需要设置该值。

***numyFine***  
类型：integer  
说明：可选项。y方向的细网格个数，单位：无。当 **meshingType** 的值为 `adaptive` 时，需要设置该值。

***numzFine***  
类型：integer  
说明：可选项。z方向的细网格个数，单位：无。当 **meshingType** 的值为 `adaptive` 时，需要设置该值。

注意：${fineLength}^{3}\times(numxFine\times numyFine\times numzFine) \le length^{3}$

#### 5.2. 网格参数配置文件

***meshType***  
类型：字符串  
说明：必选项。box的边界条件，`periodic` 或 `freeSurfaces`，默认为 `periodic`。

***length***  
类型：double precision  
说明：必选项。网格的边长，单位:nm。

***numx***  
类型：integer  
说明：必选项。x方向的细网格个数，单位：无。该值要大于等于x方向的进程数。

***numy***  
类型：integer  
说明：必选项。y方向的细网格个数，单位：无。该值要大于等于y方向的进程数。

***numz***  
类型：integer  
说明：必选项。z方向的细网格个数，单位：无。该值要大于等于z方向的进程数。


