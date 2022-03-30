%% Problem 1 {least squares approximation using SVD}
clear
x = 3;  % true slope 
A = [0: .25: 5]'; 
b = A*x + 1*randn(size(A));  %add noise 
plot(A, A*x, 'k'); % true relationship 
hold on 
plot(A,b, 'rx');

[U, S, V] = svd(A,0); % Returns the economy SVD
Xtilde = V*(inv(S))*U'*b;
plot(A, A*Xtilde, 'r');



%% Problem 2 {Ray Detection}
clear
% Data for HW 7, P2 (estimate light-ray amplitude and direction)
m = 4;
sigma = 0.5;
% theta is azimuth angle
theta =[ 3 
 10 
 80 
 150 
];
%  phi is elevation angle
phi =[ 88 
 34 
 30 
 20 
];
p = [ 1.58 
 1.50 
 2.47 
 1.10 
];


q = [];
for i = 1:length(theta)
    n = [cosd(theta(i))*cosd(phi(i)); sind(theta(i))*cosd(phi(i)); sind(phi(i))];
    q = [q n];
end


[U, S, V] = svd(q', 0); %Returns the economy SVD
w = V'*inv(S)*U'*p;

w = w/sigma;
[TH, PHI, alpha] = cart2sph(w(1),w(2),w(3));





%% Problem 3 {Closest Point}
% Data for HW 7, P3 (point closest to the given lines)
clear
n = 2;
m = 15;
P = [ -0.7930, -6.1337, 6.5463, 1.3785, -4.0749, 1.9969, 3.0305, 1.7315, -0.2904, 4.6367, 1.3919, 0.0103, 3.1012, 10.9183, 4.6092;
        2.4518, -4.4494, 2.6530, -3.6324, 3.5452, 3.8776, 0.5708, -0.0375, 1.4084, 7.7409, 3.0078, 2.4676, -1.0208, 4.2206, 6.0331;
        ];
V = [ -1.0000, -0.4409, 0.9477, -0.2011, -0.0852, 0.9502, -0.6529, -0.5113, -0.9721,0.5746, -0.5262, 0.9977, -0.5202, -0.4683, -0.6865;
        0.0058, -0.8976, 0.3193, 0.9796, -0.9964, 0.3115, -0.7575, 0.8594, -0.2347, -0.8184, -0.8504, -0.0676, 0.8540, 0.8836, -0.7272;
        ];


Q = [0, 0 ; 0, 0];
for i = 1:length(V)
    Qn = (eye(2) - (V(:,i:i)*V(:,i:i)'))^2;
    Q = Q + Qn;
end

W = [0 ; 0];
for i = 1:length(P)
    PP = ((eye(2) - (V(:,i:i)*V(:,i:i)'))^2)*(P(:,i:i));
    W = W + PP;
end

zstar = inv(Q)*W;

t=-50:0.1:50;
for i=1:m
    cor=V(:,i)*t+P(:,i)*ones(1,length(t));
    plot(cor(1,:),cor(2,:));
    hold on
end
plot(zstar(1),zstar(2),'*r')
axis([-20 20 -20 20])


