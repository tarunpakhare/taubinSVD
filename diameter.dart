import 'package:flutter/material.dart';
import 'package:scidart/numdart.dart';

class diameter {
  double calculateDiameter(List<List<Offset>> coordinates) {
    List<double> centroid = calculateCentroid(coordinates);
    List<double> X = [];
    List<double> Y = [];
    List<double> Z = [];
    double zMean = 0.0;
    int pointsCount = 0;

    for (var points in coordinates) {
      for (var point in points) {
        double deltaX = point.dx - centroid[0]; //centering data
        double deltaY = point.dy - centroid[1]; //centering data
        X.add(deltaX);
        Y.add(deltaY);
        double squaredDistance = (deltaX * deltaX) + (deltaY * deltaY); //Z = X.*X + Y.*Y;
        Z.add(squaredDistance);
        zMean += squaredDistance;
        pointsCount++;
      }
    }
  

    zMean /= pointsCount; //Zmean = mean(Z);
    print("zMeanUpdated=$zMean");
    List<double> Z0 = [];
    for (var z in Z) {
      double z0 = (z - zMean) / (2 * sqrt(zMean));
      Z0.add(z0);
    }
    print("Z0:$Z0"); //Z0 = (Z-Zmean)/(2*sqrt(Zmean));

  
    Array zo = Array(Z0);
    Array x = Array(X);
    Array y = Array(Y);
    List<Array> zxy = Array2d([zo, x, y]);
    print("ZXY:$zxy"); //ZXY = [Z0 X Y];

    Array2d zxytranspose = matrixTranspose(Array2d([zo, x, y]));
    print("Transpose:$zxytranspose");

    Array2d svdResult = svd(zxytranspose);
    print("svdResult:$svdResult"); //[U,S,V]=svd(ZXY,0);
   
    Array2d svdResultTranspose = matrixTranspose(svdResult);
    print("Transposed SVDResult:$svdResultTranspose");
    for (var i = 0; i < svdResultTranspose.length; i++) {
      for (var j = 0; j < svdResultTranspose[i].length; j++) {
        svdResultTranspose[i][j] = -svdResultTranspose[i][j];
      }
    }
    print('Changed signs: $svdResultTranspose');
    List<double> A = getcolumn(svdResultTranspose, 2); //A = V(:,3);
    print("First A:$A");
    A[0] = A[0] / (2 * sqrt(zMean)); //A(1) = A(1)/(2*sqrt(Zmean));
    print("Second A:$A");
    A.add(-(zMean) * A[0]); //A = [A ; -Zmean*A(1)];
    print("Third A:$A[0]");
    
    double disc = ((A[1] * A[1]) + (A[2] * A[2]) - (4 * A[0] * A[3]));
    print(sqrt(disc.abs()));
    print((A[0].abs() / 2));
    double radius = sqrt(disc) / A[0].abs() / 2;
    List<double> c = [];
    for (int i = 0; i < 3; i++) {
      c.add(-A[i] / A[0] / 2);
    }
    double xc = c[0] + centroid[0];
    double yc = c[1] + centroid[1];

    print("CenterX:$xc");
    print("CenterY:$yc");
    print("Radius:$radius");
    double diameter = 2 * radius;
    print("Diameter:$diameter");

    return diameter; //Par = [-(A(2:3))'/A(1)/2+centroid , sqrt(A(2)*A(2)+A(3)*A(3)-4*A(1)*A(4))/abs(A(1))/2];
  }
  //calculating the centroid
  List<double> calculateCentroid(List<List<Offset>> xy) {
    double sumX = 0.0;
    double sumY = 0.0;
    int pointCount = 0;

    for (var points in xy) {
      for (var point in points) {
        sumX += point.dx;
        sumY += point.dy;
        pointCount++;
      }
    }
    print("SumX:$sumX");
    print("SumY:$sumY");
    print("pointCount:$pointCount");
    double centerX = sumX / pointCount; 
    double centerY = sumY / pointCount; 

    print("centerX:$centerX");
    print("centerY:$centerY");

    return [centerX, centerY];
  }
  //To get column
  List<double> getcolumn(List<List<double>> matrix, int col) {
    List<double> colValues = matrix.map((row) => row[col]).toList();
    print(colValues);
    return colValues;
  }

  Array2d svd(Array2d zxy) {
    //Array2d cde = Array2d(zxy);
    SVD usv = SVD(zxy);
    Array2d u = usv.U();
    Array2d s = usv.S();
    Array2d v = usv.V();
    Array2d vtranspose = matrixTranspose(v);
    List<List<double>> ulist = u.toList();
    List<List<double>> slist = s.toList();
    List<List<double>> vlist = vtranspose.toList(growable: true);
    print("Ulist:$ulist");
    print("Slist:$slist");
    print("Vlist:$vlist");
    return vtranspose;
  }
}
