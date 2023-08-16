# KapranovDegrees


This repository contains code to compute Kapranov degrees, intersection numbers of pullbacks of psi classes on M0nbar along forgetful maps. The input in n-3 subsets of [n], S_1, ..., S_{n-3}, and i_j \in S_j for each j. We compute the Kapranov degree via
kapranovDegree([[S_1, i_1], ..., [S_{n-3}, i_{n-3}]])

For example,
kapranovDegree([[{1,2,3,4}, 1], [{1,2,3,4,5}, 2], [{3,4,5,6}, 6]])

When all of the sets of have size 4, these are called cross-ratio degrees. Then the Kapranov degree does not depend on the choice of i_j. We can compute the cross-ratio degree via
crossRatioDeg([S_1, ..., S_{n-3}])

For example,
crossRatioDeg([{1,2,3,4}, {1,2,5,6}, {3,4,5,6}])