% creating geometry for the inflow
r = 0.5;
h1 = 0.1;
l = 4;

% points
p0 = [0, 0];
p1 = [0, -r];
p2 = [0, r];
p3 = [l/4, r];
p4 = [l/8+l/4, r+h1/2];
p5 = [l/2, r + h1];
h2 = 0.5;
p6 = [l/2+l/8, r + 3*h1/2];
p7 = [13/16*l, r + 3*h1/2+h2/2];
p8 = [l, r + 3*h1/2+h2];
nnorm = sqrt( (3/16*l)*(3/16*l) + (h2/2)*(h2/2));
n1 = (3/16*l)/nnorm;
n2 = (h2/2)/nnorm;
r1 = 0.5;
t1 = n2;
t2 = -n1;
p9 = [l + t1*r1, r + 3*h1/2+h2 + t2*r1];
p10 = [13/16*l + t1*r1, r + 3*h1/2+h2/2 + t2*r1];
p11 = [l/2+l/8 + t1*r1, r + 3*h1/2 + t2*r1];
h3 = 0.00;
p13 = [l/4, -r];
p14 = [l/8+l/4, -(r+h3/2)];
p15 = [l/2, -(r + h3)];
p16 = [l/2+l/8, -(r + 3*h3/2)];
h4 = 0.1;
p17 = [13/16*l, -(r + 3*h3/2+h4/2)];
p18 = [l, -(r + 3*h3/2+h4)];
nnorm = sqrt( (3/16*l)*(3/16*l) + (h4/2)*(h4/2));
n1 = (3/16*l)/nnorm;
n2 = -(h4/2)/nnorm;
r2 = 0.3;
t1_ = -n2;
t2_ = +n1;
p19 = [l + t1_*r2, -(r + 3*h3/2+h4) + t2_*r2];
p20 = [13/16*l+ t1_*r2, -(r + 3*h3/2+h4/2)+ t2_*r2];
p21 = [l/2+l/8+ t1_*r2, -(r + 3*h3/2)+ t2_*r2];
p12 = [(l/2+l/8 + t1*r1 + l/2+l/8+ t1_*r2)/2, (r + 3*h1/2 + t2*r1 -(r + 3*h3/2)+ t2_*r2)/2];
xc = [1.75 0];

% control points for middle line
p3_ = [p3(1) 0];
p4_ = [p4(1) 0];
p5_ = [p5(1) 0];

ni = 4;
nj = 3;
coefs = zeros(2,ni,nj);

coefs(:,1,1) = p1;
coefs(:,2,1) = p13;
coefs(:,3,1) = p14;
coefs(:,4,1) = p15;

coefs(:,1,2) = p0;
coefs(:,2,2) = p3_;
coefs(:,3,2) = p4_;
coefs(:,4,2) = xc;

coefs(:,1,3) = p2;
coefs(:,2,3) = p3;
coefs(:,3,3) = p4;
coefs(:,4,3) = p5;

nurb = nrbmak(coefs,{[0 0 0 0 1 1 1 1],[0 0 0.5 1 1]});