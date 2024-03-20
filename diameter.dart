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
        double deltaX = point.dx - centroid[0];
        double deltaY = point.dy - centroid[1];
        X.add(deltaX);
        Y.add(deltaY);
        double squaredDistance = (deltaX * deltaX) + (deltaY * deltaY);
        //double z = pow(point.dx - centroid[0], 2) +
        //    pow(point.dy - centroid[1], 2).toDouble();
        Z.add(squaredDistance);
        zMean += squaredDistance;
        pointsCount++;
      }
    }
    print("Z:$Z");
    print("zMean=$zMean");

    zMean /= pointsCount; //(coordinates.length * coordinates[0].length);
    print("zMeanUpdated=$zMean");
    List<double> Z0 = [];
    for (var z in Z) {
      double z0 = (z - zMean) / (2 * sqrt(zMean));
      Z0.add(z0);
    }
    print("Z0:$Z0");

    /* List<List<double>> ZXY = [];
    for (int i = 0; i < X.length; i++) {
      ZXY.add([Z0[i].toDouble(), X[i].toDouble(), Y[i].toDouble()]);
    }
    Array2d abc = Array2d(ZXY);*/
    Array zo = Array(Z0);
    Array x = Array(X);
    Array y = Array(Y);
    List<Array> zxy = Array2d([zo, x, y]);
    print("ZXY:$zxy");

    Array2d zxytranspose = matrixTranspose(Array2d([zo, x, y]));
    print("Transpose:$zxytranspose");

    Array2d svdResult = svd(zxytranspose);
    print("svdResult:$svdResult");
    print("Till here executed successfully");
    Array2d svdResultTranspose = matrixTranspose(svdResult);
    print("Transposed SVDResult:$svdResultTranspose");
    for (var i = 0; i < svdResultTranspose.length; i++) {
      for (var j = 0; j < svdResultTranspose[i].length; j++) {
        svdResultTranspose[i][j] = -svdResultTranspose[i][j];
      }
    }
    print('Changed signs: $svdResultTranspose');
    List<double> A = getcolumn(svdResultTranspose, 2);
    print("First A:$A");
    //List<double> A = svdResult[2].toList();
    //print("A:$A");
    A[0] = A[0] / (2 * sqrt(zMean));
    print("Second A:$A");
    A.add(-(zMean) * A[0]);
    print("Third A:$A[0]");
    /*double a = A[1];
    double b = A[2];
    List<double> center = [
      -(a / (2 * A[0])) + centroid[0],
      -(b / (2 * A[0])) + centroid[1]
    ];
    print("Center:$center");*/
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

    return diameter; //aise hi rakha hai
    /* A[0] = A[0] / (2 * sqrt(zMean));
    A.add(-zMean * A[0]);

    double xCenter = -(A[1] + centroid[0] / 2 * A[0]);
    double yCenter = -(A[2] + centroid[1] / 2 * A[0]);
    double radius =
        sqrt(pow(A[1], 2)) + pow(A[2], 2) - 4 * A[0] * A[3] / (2 * A[0]);
    print(xCenter);
    print(yCenter);
    double d = 2 * radius;
    print("Diameter=$d");
    return 2 * radius;*/
  }

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
    double centerX = sumX / pointCount; //(xy.length * xy[0].length);
    double centerY = sumY / pointCount; //(xy.length * xy[0].length);

    print("centerX:$centerX");
    print("centerY:$centerY");

    return [centerX, centerY];
  }

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
