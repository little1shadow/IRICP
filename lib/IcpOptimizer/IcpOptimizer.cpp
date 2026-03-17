#include "IcpOptimizer.h"
#include <numeric>
#include "\lib\ObjLoader\ObjLoader.h"
#include<vector>
#include<chrono>
using namespace std;
using namespace Eigen;
using namespace nanoflann;

const double e = 2.7182818284;
double sigma = 2;
double sigma1=sigma;//correntropy kernel width
double sigma2=sigma;//errorentropy kernel width
double Lambda = 0.5;//the  scalar factor betwween the correntropy and   errorentropy
/*
Main constructor. Initializes the point clouds and the sparse ICP parameters.
*/
IcpOptimizer::IcpOptimizer(Matrix<double, Dynamic, 3> _firstCloud, Matrix<double, Dynamic, 3> _secondCloud, size_t _kNormals, int _nbIterations, int _nbIterationsIn, double _mu, int _nbIterShrink, double _p, IcpMethod _method, bool _verbose) : firstCloud(_firstCloud), secondCloud(_secondCloud), kNormals(_kNormals), nbIterations(_nbIterations), nbIterationsIn(_nbIterationsIn), mu(_mu), nbIterShrink(_nbIterShrink), p(_p), method(_method), verbose(_verbose)
{
  // Normal estimation
  
  
  cout << "Estimating normals for first cloud" << endl;
  //string first_path("E:/SignalProccessing/ExperimentCode/icpSparse_MEE_MCC/data_oringi/ArmadilloStand/ArmadilloStand_0.obj");
  //string first_path("E:/SignalProccessing/ExperimentCode/icpSparse_MEE_MCC/data_oringi/bunny/bunny_side2.obj");
  string first_path("E:/SignalProccessing/ExperimentCode/icpSparse_MEE_MCC/data_oringi/dragon_stand/dragonStandRight_0.obj");
  ObjectLoader myLoader;
  Matrix<double,Dynamic,3> pointCloudOne = myLoader(first_path);
  firstNormals = estimateNormals(pointCloudOne, kNormals);

  
 
  
  //firstNormals = estimateNormals(_firstCloud, kNormals);
  if (method == pointToPlane)
  {  
    cout << "Estimating normals for second cloud" << endl;
    secondNormals = estimateNormals(_secondCloud, kNormals);
    cout << "Done with normal estimation" << endl;
  }

  // Initialize the computed transformation
  computedTransfo = RigidTransfo(RotMatrix::Identity(), TransMatrix::Zero(3, 1));

  // Initialize the Lagrange multipliers to 0 for step 2.1
  lambda.resize(firstCloud.rows(), 3);
  lambda.setZero();

  // Initialize the reference distance (bounding box diagonal of cloud 1)
  Matrix<double, 1, 3> minCloudOne = firstCloud.colwise().minCoeff();
  Matrix<double, 1, 3> maxCloudOne = firstCloud.colwise().maxCoeff();
  referenceDist = (maxCloudOne - minCloudOne).norm();
  cout << "The reference distance is : " << referenceDist << endl;

  // Initialize the other parameter
  hasBeenComputed = false;
}

/*
This function is the main implementation of the algorithm where every step are made explicit.
*/
int IcpOptimizer::performIRICP()
{
  if (firstCloud.rows() == 0 || secondCloud.rows() == 0)
    return 1;
  if (method == pointToPoint)
    cout << "Beginning ICP with method Point to Point" << endl;
  else if (method == pointToPlane)
    cout << "Beginning ICP with method Point to Plane" << endl;
  // Initialize the point cloud that is going to move

  //define the groundtruth transformation
  Matrix<double,3,3> TruthRotation;
  TruthTranslation<<0,0,0;
  RigidTransfo TruthTransformation =RigidTransfo(TruthRotation,TruthTranslation);
  
  movingPC=firstCloud;
  movingNormals = firstNormals;

  

  // Beginning of the algorithm itself  nbIterations

  vector<double> ERs; //save the RMSE between movedPC and matchPC
  vector<double> Ets;

  
 
  vector<int> firstCloudGranullballPointNumber; 
  vector<double> firstCloudGranullballRadius;


  cout<<"Granullball backed run!"<<endl;

  
  ifstream input_firstCloud_GranullballPointNumber("./granuleball_info/ArmadilloStand_0_gb_info/ArmadilloStand_0_gb_info_0.78/ArmadilloStand_0_gb_0.78_pointNum.txt"); 
  ifstream input_firstCloud_GranullballRadius("./granuleball_info/ArmadilloStand_0_gb_info/ArmadilloStand_0_gb_info_0.78/ArmadilloStand_0_gb_0.78_radius.txt"); 
  //ifstream input_firstCloud_GranullballPointNumber("./granuleball_info/bunny_side2_gb_info/bunny_side2_gb_info_0.78/bunny_side2_gb_0.78_pointNum.txt"); 
  //ifstream input_firstCloud_GranullballRadius("./granuleball_info/bunny_side2_gb_info/bunny_side2_gb_info_0.78/bunny_side2_gb_0.78_radius.txt"); 
  //ifstream input_firstCloud_GranullballPointNumber("./granuleball_info/dragonStandRight_0_gb_info/dragonStandRight_0_gb_info_0.78/dragonStandRight_0_gb_0.78_pointNum.txt"); 
  //ifstream input_firstCloud_GranullballRadius("./granuleball_info/dragonStandRight_0_gb_info/dragonStandRight_0_gb_info_0.78/dragonStandRight_0_gb_0.78_radius.txt"); 


 
  if (input_firstCloud_GranullballPointNumber.is_open()) { // 如果文件成功打开
      int num;
      while (input_firstCloud_GranullballPointNumber >> num) { // 逐行读取数据
          firstCloudGranullballPointNumber.push_back(num); // 将读取的数值添加到向量中
      }
      input_firstCloud_GranullballPointNumber.close(); // 关闭文件
  }
  else {
      cout << "Unable to open source PointCloud granullball sphere number file" << endl; // 文件无法打开，输出错误信息
  }

  /*源点云粒球的半径*/
  if (input_firstCloud_GranullballRadius.is_open()) { // 如果文件成功打开
      double num;
      while (input_firstCloud_GranullballRadius >> num) { // 逐行读取数据
          firstCloudGranullballRadius.push_back(num); // 将读取的数值添加到向量中
      }
      input_firstCloud_GranullballRadius.close(); // 关闭文件
  }
  else {
      cout << "Unable to open source PointCloud granullball sphere radius file" << endl; // 文件无法打开，输出错误信息
  }

 
  vector<int> secondCloudGranullballPointNumber;
  vector<double> secondCloudGranullballRadius; 

  ifstream input_secondCloud_GranullballPointNumber("./granuleball_info/ArmadilloStand_0_gb_info/ArmadilloStand_0_gb_info_0.78/ArmadilloStand_0_rotatey_180_gb_0.78_pointNum.txt"); 
  ifstream input_secondCloud_GranullballRadius("./granuleball_info/ArmadilloStand_0_gb_info/ArmadilloStand_0_gb_info_0.78/ArmadilloStand_0_rotatey_180_gb_0.78_radius.txt"); 
  //ifstream input_secondCloud_GranullballPointNumber("./granuleball_info/bunny_side2_gb_info/bunny_side2_gb_info_0.78/bunny_side2_rotatey_180_gb_0.78_pointNum.txt"); 
  //ifstream input_secondCloud_GranullballRadius("./granuleball_info/bunny_side2_gb_info/bunny_side2_gb_info_0.78/bunny_side2_rotatey_180_gb_0.78_radius.txt"); 
  //ifstream input_secondCloud_GranullballPointNumber("./granuleball_info/dragonStandRight_0_gb_info/dragonStandRight_0_gb_info_0.78/dragonStandRight_0_rotatey_180_overlap_gb_0.78_pointNum.txt"); 
  //ifstream input_secondCloud_GranullballRadius("./granuleball_info/dragonStandRight_0_gb_info/dragonStandRight_0_gb_info_0.78/dragonStandRight_0_rotatey_180_overlap_gb_0.78_radius.txt"); 

 
  if (input_secondCloud_GranullballPointNumber.is_open()) { // 如果文件成功打开
      int num;
      while (input_secondCloud_GranullballPointNumber >> num) { // 逐行读取数据
          secondCloudGranullballPointNumber.push_back(num); // 将读取的数值添加到向量中
      }
      input_secondCloud_GranullballPointNumber.close(); // 关闭文件
  }
  else {
      cout << "Unable to open target PointCloud granullball sphere number file" << endl; // 文件无法打开，输出错误信息
  }

  /*目标点云粒球的半径*/
  if (input_secondCloud_GranullballRadius.is_open()) { // 如果文件成功打开
      double num;
      while (input_secondCloud_GranullballRadius >> num) { // 逐行读取数据
          secondCloudGranullballRadius.push_back(num); // 将读取的数值添加到向量中
      }
      input_secondCloud_GranullballRadius.close(); // 关闭文件
  }
  else {
      cout << "Unable to open target PointCloud granullball sphere radius file" << endl; // 文件无法打开，输出错误信息
  }
  cout<<"completed load granullball input file "<<endl;

  
  
  lastItertransfo = RigidTransfo(RotMatrix::Identity(), TransMatrix::Zero(3, 1));

  PointCloud source;
  PointCloud target;
  for (int iter = 0; iter < 50; iter++)
  {
    cout << "Iteration " << iter << endl<<endl;
    // 1st step : Computing correspondances
    //核宽从细到粗
    
    if((iter+1)%19==0){ 
      //sigma1=sigma1*2;
      sigma2=sigma2/2;
    }
    
    if((iter+1)%28==0){
        Lambda=1;
    }
    if((iter+1)%121==0){
      //Lambda=1;
    }
    
  
    vector<int> matchIndice = computeCorrespondances(secondCloud, movingPC);
  
    vector<pair<int, int>> PointPairWise;
    // 生成匹配对
    
    for (int i = 0; i < movingPC.rows(); i++)
    {
          pair<int, int> PointPair = {i, matchIndice[i]};
          int sourceGranullballnum = firstCloudGranullballPointNumber[i];
          int targetGranullballnum = secondCloudGranullballPointNumber[matchIndice[i]];
          double sourceGranullballRadius = firstCloudGranullballRadius[i];
          double targetGranullballRadius = secondCloudGranullballRadius[matchIndice[i]];
          double sourcedensity = (double)sourceGranullballnum/sourceGranullballRadius;
          double targetdensity = (double) targetGranullballnum/targetGranullballRadius;
          double densitydiff= abs(sourcedensity-targetdensity)/min(sourcedensity,targetdensity);
          double pointnumdiff = (double)(sourceGranullballnum - targetGranullballnum) / min(sourceGranullballnum, targetGranullballnum);
          double pointradiusdiff = (sourceGranullballRadius - targetGranullballRadius) / min(sourceGranullballRadius, targetGranullballRadius);
         
          if (densitydiff<0.2)
          {
            PointPairWise.push_back(PointPair);
          }
    }
    
    
    //查看一下每次多少正常匹配对
    cout<<"the inlier size is:"<<PointPairWise.size()<<endl;

    cout<<"the computed rotatetion is:"<<endl<<computedTransfo.first<<endl<<endl;
    cout<<"the computed translation is:"<<endl<<computedTransfo.second<<endl<<endl;
    
  
    //PointCloud matchPC = selectSubsetPC(secondCloud, matchIndice); // Selecting y
    /*生成匹配对*/
    vector<int> firstIndices(PointPairWise.size());
    vector<int> secondIndices(PointPairWise.size());
    std::transform(PointPairWise.begin(), PointPairWise.end(), firstIndices.begin(),
               [](const std::pair<int, int>& pair) {
                   return pair.first;
               });
    std::transform(PointPairWise.begin(), PointPairWise.end(), secondIndices.begin(),
               [](const std::pair<int, int>& pair) {
                   return pair.second;
               });    

    PointCloud matchPC = selectSubsetPC(secondCloud, secondIndices); 
    PointCloud firstCloudSelected =selectSubsetPC(firstCloud, firstIndices); 
    PointCloud movingPCSelected =selectSubsetPC(movingPC, firstIndices); 


    

   
    if (method == pointToPlane)
      selectedNormals = selectSubsetPC(secondNormals, matchIndice);


    // 2nd step : Computing transformation   nbIterationsIn
    RigidTransfo iterTransfo;
    for (int iterTwo = 0; iterTwo < 1; iterTwo++)
    {
     
      if (method == pointToPoint)
      {
        //iterTransfo = rigidTransformPointToPoint(firstCloudSelected, matchPC);
        auto start = std::chrono::steady_clock::now();
        iterTransfo = rigidTransformPointToPoint(movingPCSelected, matchPC);
        auto end = std::chrono::steady_clock::now();
        auto duration = end - start;

        // 输出时间差（以毫秒为单位）
        std::cout << "Time difference: " << std::chrono::duration <double, std::milli>(duration).count() << " ms" << std::endl;
      }
      
      // Updating the moving pointCloud


      //movingPC = movePointCloud(firstCloud, iterTransfo); 
      movingPC = movePointCloud(movingPC, iterTransfo);  
      movingNormals = (iterTransfo.first * movingNormals.transpose()).transpose();
      computedTransfo = compose(iterTransfo, computedTransfo);
      lastItertransfo = iterTransfo;//将最新的刚性变换更新
      updateIter(iterTransfo); // Updating the iterations measure
      

      
    }
  
    
    double tr=(TruthRotation*computedTransfo.first.transpose()).trace();
    double ER=acos((tr-1)/2);
    double Et=(TruthTranslation-computedTransfo.second).norm();
    ERs.push_back(ER);
    Ets.push_back(Et);
    source=movingPCSelected;
    target=matchPC;

  }
  vector<double> errs;
  for(int i=0;i<source.rows();i++){
    errs.push_back((source.row(i).transpose()-target.row(i).transpose()).norm());
  }
  
  
  
  
  string first_path("E:/SignalProccessing/ExperimentCode/icpSparse_MEE_MCC/data_oringi/dragon_stand/dragonStandRight_0.obj");
  ObjectLoader myLoader;
  Matrix<double,Dynamic,3> pointCloudOne = myLoader(first_path);
  movingPC=movePointCloud(pointCloudOne,computedTransfo);

  

  hasBeenComputed = true; 
  return 0;
}

/*
This function computes each closest point in refCloud for each point in queryCloud using the nanoflann kd-tree implementation. It returns the indice of the closest points of queryCloud.
*/
vector<int> IcpOptimizer::computeCorrespondances(Matrix<double, Dynamic, 3> refCloud, Matrix<double, Dynamic, 3> queryCloud)
{
  
  //Create an adapted kd tree for the point cloud
  typedef KDTreeEigenMatrixAdaptor< Matrix<double,Dynamic,3> > my_kd_tree_t;

  //Create an index
  my_kd_tree_t   mat_index(3, refCloud, 10 /* max leaf */ );
  mat_index.index->buildIndex();

  vector<int> nearestIndices;
  for(int i=0;i<queryCloud.rows();i++)
  {
    //Current query point
    Matrix<double,1,3> queryPoint = queryCloud.block(i,0,1,3);

    //Do a knn search
    const size_t num_results = 1; //We want the nearest neighbour
    vector<size_t>   ret_indexes(num_results);
    vector<double> out_dists_sqr(num_results);

    KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );

    mat_index.index->findNeighbors(resultSet, &queryPoint[0], SearchParams(10));

    //Updating the resulting index vector
    nearestIndices.push_back(ret_indexes[0]);

    if(verbose)
    {
      cout << queryPoint(0,0) << " " << queryPoint(0,1) << " " << queryPoint(0,2) << " refPoint" << endl;
      cout << refCloud(ret_indexes[0],0) << " " << refCloud(ret_indexes[0],1) << " " << refCloud(ret_indexes[0],2) << " closestPoint" << endl << endl << endl;
    }
  }
  return nearestIndices;
}



/*
Move the pointCloud according to the rigid transformation in t
*/
PointCloud IcpOptimizer::movePointCloud(PointCloud pointCloud, RigidTransfo t)
{
  return (t.first * pointCloud.transpose() + t.second.replicate(1, pointCloud.rows())).transpose();
}

/*
This function estimates the normals for the point cloud pointCloud. It makes use of the k nearest neighbour algorithm implemented in FLANN
*/
Matrix<double, Dynamic, 3> IcpOptimizer::estimateNormals(Matrix<double, Dynamic, 3> pointCloud, const size_t k)
{
  // Create an adapted kd tree for the point cloud
  typedef KDTreeEigenMatrixAdaptor<Matrix<double, Dynamic, 3>> my_kd_tree_t;

  // Create an index
  my_kd_tree_t mat_index(3, pointCloud, 10 /* max leaf */);
  mat_index.index->buildIndex();

  Matrix<double, Dynamic, 3> normals;
  normals.resize(pointCloud.rows(), 3);
  for (int i = 0; i < pointCloud.rows(); i++)
  {
    // Current point for which the normal is being computed
    Matrix<double, 1, 3> currentPoint = pointCloud.block(i, 0, 1, 3);

    // Do a knn search
    vector<size_t> ret_indexes(k);
    vector<double> out_dists_sqr(k);

    KNNResultSet<double> resultSet(k);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);

    mat_index.index->findNeighbors(resultSet, &currentPoint[0], SearchParams(10));

    // Compute the covariance matrix

    // Compute the barycentre of the nearest neighbours
    Matrix<double, 1, 3> barycentre;
    for (int j = 0; j < 3; j++)
    {
      double curVal = 0.;
      for (int neighbour = 0; neighbour < k; neighbour++)
        curVal += pointCloud(ret_indexes[neighbour], j);
      barycentre(0, j) = curVal / double(k);
    }

    // Compute the centered nearest neighbour matrix
    Matrix<double, Dynamic, 3> centeredNN;
    centeredNN.resize(k, 3);
    for (int j = 0; j < k; j++)
      centeredNN.row(j) = pointCloud.row(ret_indexes[j]) - barycentre;

    // Compute the covariance matrix
    Matrix<double, 3, 3> covariance = centeredNN.transpose() * centeredNN;

    // Computing its eigen values
    EigenSolver<Matrix<double, Dynamic, Dynamic>> eigenSolver(covariance);

    // Find the indice of the lowest eigen value
    int bestIndice = -1;
    double bestVal = DBL_MAX;
    for (int j = 0; j < 3; j++)
      if (eigenSolver.eigenvalues()(j, 0).real() < bestVal)
      {
        bestVal = eigenSolver.eigenvalues()(j, 0).real();
        bestIndice = j;
      }

    // Filling the normal
    normals.row(i) = eigenSolver.eigenvectors().block(0, bestIndice, 3, 1).normalized().transpose().real();
  }
  return normals;
}
/*
This function is the standard point to point ICP
a : moving cloud
b : reference cloud
*/
RigidTransfo IcpOptimizer::rigidTransformPointToPoint(PointCloud a, PointCloud b) const
{
  
  // the new icp based on the fusion of minimum error entropy and maxmum correntropy
  
  // compute the error of every correspondence
  vector<double> error;
  for (int i = 0; i < a.rows(); i++)
  {
    //double err = (computedTransfo.first * a.row(i).transpose() + computedTransfo.second - b.row(i).transpose()).norm();
    double err = (a.row(i).transpose() - b.row(i).transpose()).norm();//使用moving与matched进行配准
    
    if (err == 0)
    {
      error.push_back(0.0000000001);
      continue;
    }

    error.push_back(err);
  }

  // compute the weight c and w
  
  vector<double> c;
  vector<vector<double>> w;
  vector<double> w_temp; //the temp of w in each iteration
  for (int i = 0; i < a.rows(); i++)
  {
   
    for (int j = 0; j < a.rows(); j++)
    {
      //w_temp.push_back((error[i] - error[j]) * pow(e, -(error[i] - error[j]) * (error[i] - error[j]) / (2 * sigma)));
      /*使用不同核宽*/
      w_temp.push_back((-1/(sigma2*sigma2))*(error[i] - error[j]) * pow(e, -(error[i] - error[j]) * (error[i] - error[j]) / (2 * sigma2*sigma2)));
    }
    w.push_back(w_temp);
    
    w_temp.clear();
    /*不同核宽*/
    c.push_back(-pow(e, -error[i] * error[i] / (2 * sigma1*sigma1))/(sigma1*sigma1));
    
    
  }
  
  //compute W;
  double W=0;
  for (int i = 0; i < a.rows(); i++)
  {
    for (int j = 0; j < b.rows(); j++)
    {
      W = W + Lambda * c[i] + (1-Lambda) * w[i][j] * (1 / error[i] - 1 / error[j]);
      w[i][j] = (1-Lambda) * w[i][j]; // compute new w
    }
    c[i]=Lambda * c[i];//compute new c
   
  }
  W=-W;//W好像要取个负数
  // compute S
  Matrix<double, 3, 1> Sx = Matrix<double, 3, 1>::Zero(3, 1);
  Matrix<double, 3, 1> Sy = Matrix<double, 3, 1>::Zero(3, 1);

  for (int i = 0; i < a.rows(); i++)
  {
    for (int j = 0; j < b.rows(); j++)
    {
      Sx = Sx + c[i] * a.row(i).transpose() +w[i][j] * (a.row(i).transpose() / error[i] - a.row(j).transpose() / error[j]);
      Sy = Sy + w[i][j] * (b.row(j).transpose() / error[j] - b.row(i).transpose() / error[i]) - c[i] * b.row(i).transpose();
    }
  }
 
  // compute p,q
  PointCloud p = a;
  PointCloud q = b;
  for (int i = 0; i < a.rows(); i++)
  {
    p.row(i) = a.row(i) + Sx.transpose() / W;
    q.row(i) = b.row(i) - Sy.transpose() / W;
  }

  //compute u and v
  vector<double> u;
  vector<vector<double>> v;
  
  for (int i = 0; i < a.rows(); i++)
  {
    vector<double> v_temp;
    for (int j = 0; j < b.rows(); j++)
    { 
      /*逐步迭代的方式*/
      v_temp.push_back(-pow(e, -pow(((lastItertransfo.first * p.row(i).transpose() - q.row(i).transpose()).norm() - (lastItertransfo.first * p.row(j).transpose() - q.row(j).transpose()).norm()), 2) / (2 * sigma2*sigma2)));
      
      /*不同核宽*/
      //v_temp.push_back(-pow(e, -pow(((computedTransfo.first * p.row(i).transpose() - q.row(i).transpose()).norm() - (computedTransfo.first * p.row(j).transpose() - q.row(j).transpose()).norm()), 2) / (2 * sigma2*sigma2)));
    } 
    v.push_back(v_temp);
    /*使用逐步迭代的方式*/
    
    u.push_back(-pow(e, -pow((lastItertransfo.first * p.row(i).transpose() - q.row(i).transpose()).norm(),2) / (2 * sigma1*sigma1)));
    
    /*使用不同的核宽*/
    //u.push_back(-pow(e, -pow((computedTransfo.first * p.row(i).transpose() - q.row(i).transpose()).norm(),2) / (2 * sigma1*sigma1)));
    v_temp.clear();
  }

  //compute H1
  Matrix<double, 3, 3> H1 = Matrix<double, 3, 3>::Zero(3, 3);
  for (int i = 0; i < a.rows(); i++)
  {
    H1 = H1 + u[i] * p.row(i).transpose() * q.row(i);
  }
  /*不同核宽*/
  H1 = Lambda * H1 / (a.rows() * sigma1*sigma1);
  
  //compute vi and vj
  
  vector<double> vi;
  vector<double> vj;

  for (int i = 0; i < a.rows(); i++)
  {
    vi.push_back(accumulate(v[i].begin(), v[i].end(), 0.0));
  }
 
  for (int j = 0; j < a.rows(); j++)
  {
    double temp = 0;
    for (int i = 0; i < a.rows(); i++)
    {
      temp = temp + v[i][j];
    }
    vj.push_back(temp);
  }
  
  //compute m
  vector<double> m;
  for(int i = 0;i<a.rows();i++){
    m.push_back(vi[i]+vj[i]);
  }
  

  //compute H2
  Matrix<double, 3, 3> H2 = Matrix<double, 3, 3>::Zero(3, 3);
  for (int i = 0; i < a.rows(); i++)
  {
    H2 = H2 + m[i] * p.row(i).transpose() * q.row(i);
  }
  /*不同核宽*/
  H2 =(1 - Lambda) * H2 / (a.rows() * a.rows() * sigma2*sigma2);

  //compute S
  Matrix<double, 3, 3> S = Matrix<double, 3, 3>::Zero(3, 3);
  S = H1+ H2;

  //compute the svd of S
  JacobiSVD<Matrix<double, 3, 3>> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);

  // compyrte D
  Matrix<double, 3, 3> D = Matrix<double, 3, 3>::Identity(3, 3);
  if (S.determinant() > 0)
  {
    D = -1 * D;
    D.row(2) = -D.row(2);
  }
  else
  {
    D = -D;
  }

  // compute rotation
  RotMatrix rotation = svd.matrixV() * D * svd.matrixU().transpose();

  Matrix<double, 3, 1> trans = Matrix<double, 3, 1>::Zero(3, 1);
  for (int i = 0; i < a.rows(); i++)
  {
    for (int j = 0; j < a.rows(); j++)
    {
      trans = trans +c[i]*(rotation*a.row(i).transpose()-b.row(i).transpose())+ w[i][j] * ((rotation * a.row(i).transpose()-b.row(i).transpose()) / error[i] - (rotation * a.row(j).transpose()-b.row(j).transpose()) / error[j]);
    }
  }

  TransMatrix translation = trans / W;
  
  
  // Outputing the transformation
  
  if (verbose)
  {
    cout << endl
         << endl
         << "Rotation Matrix : " << endl
         << rotation << endl;
    cout << "Translation Matrix : " << endl
         << translation << endl
         << endl
         << endl;
  }
  

  return RigidTransfo(rotation, translation);
}

/*
This function is the standard point to plane ICP
a : moving cloud
b : reference cloud
n : normal to the reference cloud b
*/
RigidTransfo IcpOptimizer::rigidTransformPointToPlane(PointCloud a, PointCloud b, Matrix<double, Dynamic, 3> n) const
{
  // Initialize linear system
  Matrix<double, 6, 6> leftMember = Matrix<double, 6, 6>::Zero(6, 6);
  Matrix<double, 6, 1> rightMember = Matrix<double, 6, 1>::Zero(6, 1);

  // PointCloud c = PointCloud::Zero(a.rows(),3);
  for (int i = 0; i < a.rows(); i++)
  {
    // Computing c = a x n
    Matrix<double, 1, 3> c = a.row(i).cross(n.row(i));

    // Updating left member
    leftMember.block(0, 0, 3, 3) += c.transpose() * c;               // Top-left block
    leftMember.block(3, 3, 3, 3) += n.row(i).transpose() * n.row(i); // Bottom-right block
    leftMember.block(0, 3, 3, 3) +=
        n.row(i).replicate(3, 1).cwiseProduct(c.transpose().replicate(1, 3)); // Top-right block
    leftMember.block(3, 0, 3, 3) +=
        n.row(i).transpose().replicate(1, 3).cwiseProduct(c.replicate(3, 1)); // Bottom-left block

    // Updating right member
    double factor = (a.row(i) - b.row(i)) * n.row(i).transpose();
    rightMember.block(0, 0, 3, 1) -= factor * c.transpose();        // Top 3 elements
    rightMember.block(3, 0, 3, 1) -= factor * n.row(i).transpose(); // Bottom 3 elements
  }

  // Solving linear system
  LDLT<Matrix<double, 6, 6>> ldlt(leftMember);
  Matrix<double, 6, 1> solution = ldlt.solve(rightMember);

  // Expressing the resulting transformation
  RotMatrix rotation = (AngleAxisd(AngleAxisd::Scalar(solution(0, 0)), Vector3d::UnitX()) * AngleAxisd(AngleAxisd::Scalar(solution(1, 0)), Vector3d::UnitY()) * AngleAxisd(AngleAxisd::Scalar(solution(2, 0)), Vector3d::UnitZ())).matrix();
  TransMatrix translation = solution.block(3, 0, 3, 1);

  return RigidTransfo(rotation, translation);
}

/*
This function implements the shrink operator which optimizes the function
f(z) = ||z||_2^p + mu/2*||z-h||_2^2
*/
TransMatrix IcpOptimizer::shrink(TransMatrix h) const
{
  double alpha_a = pow((2. / mu) * (1. - p), 1. / (2. - p));
  double hTilde = alpha_a + (p / mu) * pow(alpha_a, p - 1);
  double hNorm = h.norm();
  if (hNorm <= hTilde)
    return 0 * h;
  double beta = ((alpha_a) / hNorm + 1.) / 2.;
  for (int i = 0; i < nbIterShrink; i++)
    beta = 1 - (p / mu) * pow(hNorm, p - 2.) * pow(beta, p - 1);
  return beta * h;
}

/*
Computing composition of tNew by tOld (tNew o tOld)
*/
RigidTransfo IcpOptimizer::compose(RigidTransfo tNew, RigidTransfo tOld) const
{
  return RigidTransfo(tNew.first * tOld.first, tNew.first * tOld.second + tNew.second);
  //return tNew;
}

/*
Selects the subset of rows whose index is in indice in the Point Cloud p
*/
PointCloud IcpOptimizer::selectSubsetPC(PointCloud p, vector<int> indice) const
{
  PointCloud selection = PointCloud::Zero(indice.size(), 3);
  for (int i = 0; i < indice.size(); i++)
    selection.row(i) = p.row(indice[i]);
  return selection;
}

/*
Updates the iterations measure by estimating the amplitude of rigid motion t
*/
void IcpOptimizer::updateIter(RigidTransfo t)
{
  Matrix<double, 4, 4> id = Matrix<double, 4, 4>::Identity();
  Matrix<double, 4, 4> curT = Matrix<double, 4, 4>::Identity();
  curT.block(0, 0, 3, 3) = t.first;
  curT.block(0, 3, 3, 1) = t.second / referenceDist;
  Matrix<double, 4, 4> diff = curT - id;                   // Difference between t and identity ，计算当前变换与初始变换的差异，初始变换为单位阵
  iterations.push_back((diff * diff.transpose()).trace()); // Returning matrix norm
}

/*
Save iterations to file
*/
void IcpOptimizer::saveIter(string pathToFile)
{
  ofstream txtStream(pathToFile.c_str());
  for (int i = 0; i < iterations.size(); i++)
    txtStream << iterations[i] << endl;
  txtStream.close();
}

/*
Just a getter to the normals of the first cloud (moving cloud)
*/
Matrix<double, Dynamic, 3> IcpOptimizer::getFirstNormals() const
{
  return firstNormals;
}

/*
Return a copy of the first point cloud which has been moved by the computed
rigid motion.
If the rigid motion has not been computed it returns just the original first point cloud.
*/
PointCloud IcpOptimizer::getMovedNormals() const
{
  if (hasBeenComputed)
    return firstNormals;
  else
  {
    cout << "Warning ! The transformation has not been computed ! Please use the method \
    performSparceICP() before retrieving the moved normals."
         << endl;
    return movingNormals;
  }
}

/*
Return a copy of the first point cloud which has been moved by the computed
rigid motion.
If the rigid motion has not been computed it returns just the original first point cloud.
*/
PointCloud IcpOptimizer::getMovedPointCloud() const
{
  if (hasBeenComputed)
    return movingPC;
  else
  {
    cout << "Warning ! The transformation has not been computed ! Please use the method \
    performSparceICP() before retrieving the moved point cloud."
         << endl;
    return firstCloud;
  }
}

/*
Return the computed transformation.
If it has not been computed, just returns the identity.
*/
RigidTransfo IcpOptimizer::getComputedTransfo() const
{
  if (!hasBeenComputed)
    cout << "Warning ! The transformation has not been computed ! Please use the method \
      performSparceICP() before retrieving the rigid motion."
         << endl;
  return computedTransfo;
}

/*
Returns the reference distance which is the length of the great diagonal of the first
point cloud's bounding box.
*/
double IcpOptimizer::getReferenceDist() const
{
  return referenceDist;
}

double IcpOptimizer::computeRMSE(PointCloud a, PointCloud b)
{
  double rmse = 0;
  for (int i = 0; i < a.rows(); i++)
  {
    rmse = rmse + pow((a.row(i) - b.row(i)).norm(), 2);
  }
  return rmse / a.rows();
}
// computeObjectiveFunctionValue
double IcpOptimizer::computeObjectiveFunctionValue(PointCloud a, PointCloud b)
{
  double ofv = 0;
  for (int i = 0; i < a.rows(); i++)
  {
    // correntropy
    ofv = ofv + pow(e, -pow((a.row(i) - b.row(i)).norm(), 2) / (2 * sigma));
    // p_norm
    // ofv = ofv+-pow((a.row(i)-b.row(i)).norm(),p);
  }
  return ofv / a.rows();
}

vector<double> IcpOptimizer::computedTransformationError(RigidTransfo estimation, RigidTransfo truth)
{
  vector<double> errors;
  errors.push_back(acos(((truth.first * estimation.first.transpose()).trace() - 1) / 2));
  errors.push_back((estimation.second - truth.second).lpNorm<1>());
  return errors;
}
void IcpOptimizer::saveMetrics(std::string pathToFile, vector<double> metrics)
{
  ofstream txtStream(pathToFile.c_str());
  for (int i = 0; i < metrics.size(); i++)
    txtStream << metrics[i] << endl;
  txtStream.close();
}
/*
this function compute the fpfh descriptor of each point
*/
Eigen::Matrix<double,Eigen::Dynamic,10*4> IcpOptimizer::computedFpfhDescriptor(PointCloud a,int k){
   //计算点云的FPFH描述符
   Eigen::Matrix<double, Eigen::Dynamic, 3> refCloud = a;
  
  // 计算点云的每个点的法向量
  Eigen::Matrix<double, Eigen::Dynamic, 3> sourceNormals = estimateNormals(refCloud, 10);

  // Create an adapted kd tree for the point cloud
  typedef KDTreeEigenMatrixAdaptor<Matrix<double, Dynamic, 3>> my_kd_tree_t; // 不写第三个参数，默认用欧式距离

  // Create an index

  my_kd_tree_t source_mat_index(3, refCloud, 11 );
  source_mat_index.index->buildIndex();
  

  vector<vector<size_t>> sourceNighbours;

  
  // 查找source中每个点的邻居
  //  对于每个点
  for (int i = 0; i < refCloud.rows(); ++i)
  {
    Matrix<double,1,3> queryPoint = refCloud.block(i,0,1,3);

    //Do a knn search
    const size_t num_results = k; //We want the 11 nearest neighbours
    vector<size_t>   ret_indexes(num_results);
    vector<double> out_dists_sqr(num_results);

    KNNResultSet<double> resultSet(num_results);
    resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );

    source_mat_index.index->findNeighbors(resultSet, &queryPoint[0], SearchParams(11));

    //保存邻近点下标，但是不要自己。
    vector<size_t> sourcePointNighbourIndices;
    for (int j = 1; j <  ret_indexes.size(); j++)
    {
      sourcePointNighbourIndices.push_back(ret_indexes[j]);
    }
    sourceNighbours.push_back(sourcePointNighbourIndices);
    sourcePointNighbourIndices.clear();
  }
  Eigen::VectorXd p(3); // 将点云坐标转换为向量
  Eigen::VectorXd q(3);
  vector<vector<double>> sourceDescriptor;
  // 计算原点云的描述符SPFH
  for (int i = 0; i < refCloud.rows(); i++)
  {
    vector<double> pointDescriptor;  
    for (int j = 0; j < sourceNighbours[i].size(); j++)
    { 
      p = refCloud.row(i);
      q = refCloud.row(sourceNighbours[i][j]);
    
      double distance = (p - q).norm();
      Eigen::Vector3d pq = (p - q).head<3>() / distance;
      
      Eigen::Vector3d u = sourceNormals.row(i);
      
      Eigen::Vector3d v = pq.cross(u);
      
      Eigen::Vector3d w = u.cross(v);
    
      double alpha = v.dot(sourceNormals.row(j));
      
      double phi = u.dot(pq);
      ;
      double theta = atan(w.dot(sourceNormals.row(j)) / w.dot(u));
     
      vector<double> jthdecriptor{alpha, phi, theta, distance};
      
      for (int s = 0; s < jthdecriptor.size(); s++)
      {
        pointDescriptor.push_back(jthdecriptor[s]);
      }
      
    }
    sourceDescriptor.push_back(pointDescriptor);
    pointDescriptor.clear();
  }
  // 计算源点云的FPFH
  vector<vector<double>> sourceFpfhDescriptor;
  for (int i = 0; i < sourceDescriptor.size(); i++)
  {
    vector<double> FPFH = sourceDescriptor[i];
    for (int j = 0; j < sourceNighbours[i].size(); j++)
    {
      p = refCloud.row(i);
      q = refCloud.row(sourceNighbours[i][j]);
      double distance = (p - q).norm();
      double w = 1 / (distance);
      for (int k = 0; k < FPFH.size(); k++)
      {
        FPFH[k] += w * sourceDescriptor[sourceNighbours[i][j]][k];
      }
    }
    sourceFpfhDescriptor.push_back(FPFH);
    FPFH.clear();
  }
  // 将描述符转换为 矩阵格式
  
  
  Eigen::Matrix<double, Eigen::Dynamic, 10 * 4> sourcePointCloudFpfhDescriptor; // 源点云描述符矩阵格式  /*这里还有点问题*/
  sourcePointCloudFpfhDescriptor.resize(sourceFpfhDescriptor.size(), 10 * 4);
  for (int i = 0; i < sourceFpfhDescriptor.size(); i++)
  { 
    for (int j = 0; j < sourceFpfhDescriptor[i].size(); j++)
    {
      sourcePointCloudFpfhDescriptor(i,j) = sourceFpfhDescriptor[i][j];
    }
  }

  return sourcePointCloudFpfhDescriptor;

}
