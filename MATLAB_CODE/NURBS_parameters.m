function [NURBS,Boundary,name] = NURBS_parameters(nro)

if nro == 1

    degree = 3;
    controlPoints = [ 0 -0.2 -0.5 -0.7 -1.7 -4 -4 -3.75 0 4 4 4  2.5 0;...
                     -2 -1.8 -1.3 -1.4 0.2 -3.5  0  3.5 4 4 0 -4 -0.1 -2];
    knots = [0 0 0 0 0.25 0.3 0.4 0.5 0.6 0.75 0.75 0.75 0.8 0.9 1 1 1 1];
    weights = [1 0.707 1 0.707 1 0.5 1 0.707 1 0.707 1 1 2 1];

    Boundary = [-6  6 6 -6;...
                -6 -6 6 6];

    name = "Mobi-Dick";

elseif nro == 2

    degree = 2;
    controlPoints = [ 0 -1 -1 -1 0 1 1 1 0;...
                     -1 -1 0 1 1 1 0 -1 -1];
    knots = [0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1];
    weights = [1 0.707 1 0.707 1 0.707 1 0.707 1];

    Boundary = [-5 5 5 -5;...
                -5 -5 5 5];

    name = "Circumference";

elseif nro == 3

    degree = 2;
    controlPoints = [ 0 -4 -4 -4 0 4 4 4 0;...
                     -0.9 -4 0 4 -0.3 4 0 -4 -0.9];
    knots = [0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1];
    weights = [1 0.707 1 0.707 1 0.707 1 0.707 1];

    Boundary = [-5 5 5 -5;...
                -5 -5 5 5];

    name = "DoubleCirumference";

elseif nro == 4

    degree = 3;
    controlPoints = [-2 -10 -20  -20 -10 0 10 20 20 10 0 -2;...
                      0 -2   -1    1  2  3  2  1 -1 -2 1  0];
    knots = [0 0 0 0 0.267 0.333 0.4 0.467 0.5333 0.6 0.667 0.733 1 1 1 1];
    weights = [1 0.707 1 0.707 1 0.707 1 0.707 1 0.707 1 1];

    Boundary = [-25 25 25 -25;...
                -25 -25 25 25];

    name = "FlatShape";

elseif nro == 5

    degree = 3;
    controlPoints = [ 0 -5 -2 0 2 5  0;...
                     -3.5  1.5  4.5 2.5 4.5 1.5 -3.5];
    knots = [0 0 0 0 0.5 0.5 0.5 1 1 1 1];
    weights = [1 2 2 1 2 2 1];
    controlPoints = [10/12 0; 0 10/12]*controlPoints;

    Boundary=[-5 5 5 -5;...
              -5 -5 5 5];

    name = "Heart";

elseif nro == 6

    degree = 2;
    controlPoints = [ 0 -1 -1 -1 0 1 1 1 0;...
                     -1 -1 0 1 1 1 0 -1 -1]*3;
    knots = [0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1];
    weights = [1 0.707 1 0.707 1 0.707 1 0.707 1];

    Boundary = [0 5 5 0;...
                0 0 5 5];

    name = "QuarterCircumference";

end

    NURBS.degree = degree;
    NURBS.knots = knots;
    NURBS.controlPoints = controlPoints;
    NURBS.weights = weights;

end
 
