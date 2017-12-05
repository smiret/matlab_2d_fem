clear all
close all
%Define the Parameters
kk = 0.04;
z = 500;
q = -25;
Tbar = 100;
rho = 2;
cp = 2;
ti = 0;
tf = 100;
steps = 20000;
deltat = tf/steps;

Nr = 10;
Ntheta = 20;
rshort = 0.1;
rlong = 0.25;
%Get meshing data from mesher
[conn, nodecart,nodepolar]=mesh2(rlong,rshort,Nr,Ntheta);
Ne = length(conn);
NP=length(nodecart(:,1)); % number of points

K=sparse(zeros(NP,NP)); 
M=sparse(zeros(NP,NP)); 
M2=sparse(zeros(NP,NP));
R=zeros(NP,1);

%% Take in node input
for e = 1:1:Ne
x1 = [nodecart(conn(e,1),1) nodecart(conn(e,2),1) nodecart(conn(e,3),1) nodecart(conn(e,4),1)]; 
x2 = [nodecart(conn(e,1),2) nodecart(conn(e,2),2) nodecart(conn(e,3),2) nodecart(conn(e,4),2)];

%Gauss Vector
w5 = [0.568888888888889 0.478628670499366 0.478628670499366 0.236926885056189 0.236926885056189];
gauss5 = [0.000000000000000 0.538469310105683 -0.538469310105683 0.906179845938664 -0.906179845938664];



%% Define the components for the master element
%Master Element Shape Functions
phat1 =@(z1,z2) (1/4).*(1-z1).*(1-z2); phat2 =@(z1,z2) (1/4).*(1+z1).*(1-z2);
phat3 =@(z1,z2) (1/4).*(1+z1).*(1+z2); phat4 =@(z1,z2) (1/4).*(1-z1).*(1+z2);
 

%Master Element Shape Functions Derivatives
dphat1dz1 =@(z1,z2) -(1/4)*(1-z2); dphat2dz1 =@(z1,z2) (1/4)*(1-z2);
dphat3dz1 =@(z1,z2) (1/4)*(1+z2); dphat4dz1 =@(z1,z2) -(1/4)*(1+z2);


dphat1dz2 =@(z1,z2) -(1/4)*(1-z1); dphat2dz2 =@(z1,z2) -(1/4)*(1+z1);
dphat3dz2 =@(z1,z2) (1/4)*(1+z1); dphat4dz2 =@(z1,z2) (1/4)*(1-z1);

%Mapping to Master Domain

xz1 =@(z1,z2) x1(1)*phat1 + x1(2)*phat2 + x1(3)*phat3 + x1(4)*phat4;
xz2 =@(z1,z2) x2(1)*phat1 + x2(2)*phat2 + x2(3)*phat3 + x2(4)*phat4;


%% Elemental K matrix


%First Component
 for A = 1:4
     for B = 1:4
         count = 1;
         inter = zeros(25,1);
        for ii = 1:5
            for jj = 1:5
         
                F11 = x1(1)*dphat1dz1(gauss5(ii),gauss5(jj)) + x1(2)*dphat2dz1(gauss5(ii),gauss5(jj)) + x1(3)*dphat3dz1(gauss5(ii),gauss5(jj)) + x1(4)*dphat4dz1(gauss5(ii),gauss5(jj));
                F12 = x1(1)*dphat1dz2(gauss5(ii),gauss5(jj)) + x1(2)*dphat2dz2(gauss5(ii),gauss5(jj)) + x1(3)*dphat3dz2(gauss5(ii),gauss5(jj)) + x1(4)*dphat4dz2(gauss5(ii),gauss5(jj));
                F21 = x2(1)*dphat1dz1(gauss5(ii),gauss5(jj)) + x2(2)*dphat2dz1(gauss5(ii),gauss5(jj)) + x2(3)*dphat3dz1(gauss5(ii),gauss5(jj)) + x2(4)*dphat4dz1(gauss5(ii),gauss5(jj));
                F22 = x2(1)*dphat1dz2(gauss5(ii),gauss5(jj)) + x2(2)*dphat2dz2(gauss5(ii),gauss5(jj)) + x2(3)*dphat3dz2(gauss5(ii),gauss5(jj)) + x2(4)*dphat4dz2(gauss5(ii),gauss5(jj));
                
                F = [F11 F12;
                     F21 F22];
                
                 dphatdz1vec = [dphat1dz1(gauss5(ii),gauss5(jj)) dphat2dz1(gauss5(ii),gauss5(jj)) dphat3dz1(gauss5(ii),gauss5(jj)) dphat4dz1(gauss5(ii),gauss5(jj))];
                 dphatdz2vec = [dphat1dz2(gauss5(ii),gauss5(jj)) dphat2dz2(gauss5(ii),gauss5(jj)) dphat3dz2(gauss5(ii),gauss5(jj)) dphat4dz2(gauss5(ii),gauss5(jj))];
                 
                 
                 inter1 = (transpose(inv(F))*[dphatdz1vec(A); dphatdz2vec(A)]);
                 inter2 = (transpose(inv(F))*[dphatdz1vec(B); dphatdz2vec(B)]);
                 
                 
                 inter(count) = w5(ii)*w5(jj)*dot(inter1,inter2)*kk*(det(F));
                 
         fK1(A,B) = sum(inter);
           
           count = count+1;
           
            end
        end
     end
 end

 Kstore1(e) = {fK1};
 
K(conn(e,:),conn(e,:)) = K(conn(e,:),conn(e,:)) + fK1;




end

%Boundary Condition for K
for e=1:Nr:Ne
% 
% 
 K(conn(e,1),:) = 0;
 
 K(conn(e,4),:) = 0;
 
 K(conn(e,1),conn(e,1)) = 1;
 
 K(conn(e,4),conn(e,4)) = 1;
 
 
 end


%% Construction of the R vector

%First Component in Elemental Form
for e = 1:1:Ne
x1 = [nodecart(conn(e,1),1) nodecart(conn(e,2),1) nodecart(conn(e,3),1) nodecart(conn(e,4),1)]; 
x2 = [nodecart(conn(e,1),2) nodecart(conn(e,2),2) nodecart(conn(e,3),2) nodecart(conn(e,4),2)];

%Master Element Shape Functions
phat1 =@(z1,z2) (1/4).*(1-z1).*(1-z2); phat2 =@(z1,z2) (1/4).*(1+z1).*(1-z2);
phat3 =@(z1,z2) (1/4).*(1+z1).*(1+z2); phat4 =@(z1,z2) (1/4).*(1-z1).*(1+z2);
 

%Master Element Shape Functions Derivatives
dphat1dz1 =@(z1,z2) -(1/4)*(1-z2); dphat2dz1 =@(z1,z2) (1/4)*(1-z2);
dphat3dz1 =@(z1,z2) (1/4)*(1+z2); dphat4dz1 =@(z1,z2) -(1/4)*(1+z2);


dphat1dz2 =@(z1,z2) -(1/4)*(1-z1); dphat2dz2 =@(z1,z2) -(1/4)*(1+z1);
dphat3dz2 =@(z1,z2) (1/4)*(1+z1); dphat4dz2 =@(z1,z2) (1/4)*(1-z1);

for A = 1:4
    count = 1;
    inter = zeros(25,1);
    for ii = 1:5
            for jj = 1:5
                
                F11 = x1(1)*dphat1dz1(gauss5(ii),gauss5(jj)) + x1(2)*dphat2dz1(gauss5(ii),gauss5(jj)) + x1(3)*dphat3dz1(gauss5(ii),gauss5(jj)) + x1(4)*dphat4dz1(gauss5(ii),gauss5(jj));
                F12 = x1(1)*dphat1dz2(gauss5(ii),gauss5(jj)) + x1(2)*dphat2dz2(gauss5(ii),gauss5(jj)) + x1(3)*dphat3dz2(gauss5(ii),gauss5(jj)) + x1(4)*dphat4dz2(gauss5(ii),gauss5(jj));
                F21 = x2(1)*dphat1dz1(gauss5(ii),gauss5(jj)) + x2(2)*dphat2dz1(gauss5(ii),gauss5(jj)) + x2(3)*dphat3dz1(gauss5(ii),gauss5(jj)) + x2(4)*dphat4dz1(gauss5(ii),gauss5(jj));
                F22 = x2(1)*dphat1dz2(gauss5(ii),gauss5(jj)) + x2(2)*dphat2dz2(gauss5(ii),gauss5(jj)) + x2(3)*dphat3dz2(gauss5(ii),gauss5(jj)) + x2(4)*dphat4dz2(gauss5(ii),gauss5(jj));
                
                F = [F11 F12;
                     F21 F22];
                
                 phatvec = [phat1(gauss5(ii),gauss5(jj)) phat2(gauss5(ii),gauss5(jj)) phat3(gauss5(ii),gauss5(jj)) phat4(gauss5(ii),gauss5(jj))];
                
                 inter1 = phatvec(A)*z;
                 
         inter(count) = w5(ii)*w5(jj)*(inter1*det(F));
         
         fR11(A) = sum(inter);
         
         count = count+1;
           
        end
    end
end

Rstore1(e) = {fR11}; 

fR111(:,1) = fR11;
% fR112(:,1) = fR12;
% fR113(:,1) = fR13;
% fR114(:,1) = fR14;

R(conn(e,:)) = R(conn(e,:)) + fR111;


end


%Second, flux component of R
for e = Nr:Nr:Ne
x1 = [nodecart(conn(e,1),1) nodecart(conn(e,2),1) nodecart(conn(e,3),1) nodecart(conn(e,4),1)]; 
x2 = [nodecart(conn(e,1),2) nodecart(conn(e,2),2) nodecart(conn(e,3),2) nodecart(conn(e,4),2)];

%Master Element Shape Functions
phat1 =@(z1,z2) (1/4).*(1-z1).*(1-z2); phat2 =@(z1,z2) (1/4).*(1+z1).*(1-z2);
phat3 =@(z1,z2) (1/4).*(1+z1).*(1+z2); phat4 =@(z1,z2) (1/4).*(1-z1).*(1+z2);
 

%Master Element Shape Functions Derivatives
dphat1dz1 =@(z1,z2) -(1/4)*(1-z2); dphat2dz1 =@(z1,z2) (1/4)*(1-z2);
dphat3dz1 =@(z1,z2) (1/4)*(1+z2); dphat4dz1 =@(z1,z2) -(1/4)*(1+z2);


dphat1dz2 =@(z1,z2) -(1/4)*(1-z1); dphat2dz2 =@(z1,z2) -(1/4)*(1+z1);
dphat3dz2 =@(z1,z2) (1/4)*(1+z1); dphat4dz2 =@(z1,z2) (1/4)*(1-z1);


for A = 1:4
    count = 1;
    inter = zeros(5,1);
    for jj = 1:5

         
                zz1 = 1;
                N = [zz1; 0];
                
                F11 = [x1(1)*dphat1dz1(zz1,gauss5(jj)) + x1(2)*dphat2dz1(zz1,gauss5(jj)) + x1(3)*dphat3dz1(zz1,gauss5(jj)) + x1(4)*dphat4dz1(zz1,gauss5(jj))];
                F12 = [x1(1)*dphat1dz2(zz1,gauss5(jj)) + x1(2)*dphat2dz2(zz1,gauss5(jj)) + x1(3)*dphat3dz2(zz1,gauss5(jj)) + x1(4)*dphat4dz2(zz1,gauss5(jj))];
                F21 = [x2(1)*dphat1dz1(zz1,gauss5(jj)) + x2(2)*dphat2dz1(zz1,gauss5(jj)) + x2(3)*dphat3dz1(zz1,gauss5(jj)) + x2(4)*dphat4dz1(zz1,gauss5(jj))];
                F22 = [x2(1)*dphat1dz2(zz1,gauss5(jj)) + x2(2)*dphat2dz2(zz1,gauss5(jj)) + x2(3)*dphat3dz2(zz1,gauss5(jj)) + x2(4)*dphat4dz2(zz1,gauss5(jj))];
                
                F = [F11 F12;
                     F21 F22];
                
                 phatvec = [phat1(zz1,gauss5(jj)) phat2(zz1,gauss5(jj)) phat3(zz1,gauss5(jj)) phat4(zz1,gauss5(jj))];
                 
                 inter1 = phatvec(A)*q;
                 inter2 = (dot(N,inv(F)*transpose(inv(F))*N))^0.5;
                 
         inter(count) = w5(jj)*(inter1*det(F)*inter2);
         
         fR2(A) = sum(inter);
        
         count = count+1;
    end
end

Rstore2(e) = {fR2};

fR22(:,1) = fR2;

R(conn(e,:)) = R(conn(e,:)) + fR22;

end




%Boundary Condition for R

for e = 1:Nr:Ne


R(conn(e,1)) = Tbar;
R(conn(e,4)) = Tbar;

end

%Third Component in Elemental Form
%Surface at zeta1 = -1



%% Elemental M matrix

for e= 1:1:Ne
x1 = [nodecart(conn(e,1),1) nodecart(conn(e,2),1) nodecart(conn(e,3),1) nodecart(conn(e,4),1)]; 
x2 = [nodecart(conn(e,1),2) nodecart(conn(e,2),2) nodecart(conn(e,3),2) nodecart(conn(e,4),2)];

%Master Element Shape Functions
phat1 =@(z1,z2) (1/4).*(1-z1).*(1-z2); phat2 =@(z1,z2) (1/4).*(1+z1).*(1-z2);
phat3 =@(z1,z2) (1/4).*(1+z1).*(1+z2); phat4 =@(z1,z2) (1/4).*(1-z1).*(1+z2);
 

%Master Element Shape Functions Derivatives
dphat1dz1 =@(z1,z2) -(1/4)*(1-z2); dphat2dz1 =@(z1,z2) (1/4)*(1-z2);
dphat3dz1 =@(z1,z2) (1/4)*(1+z2); dphat4dz1 =@(z1,z2) -(1/4)*(1+z2);


dphat1dz2 =@(z1,z2) -(1/4)*(1-z1); dphat2dz2 =@(z1,z2) -(1/4)*(1+z1);
dphat3dz2 =@(z1,z2) (1/4)*(1+z1); dphat4dz2 =@(z1,z2) (1/4)*(1-z1);
    
    
    
 for A = 1:4
     for B = 1:4
         count = 1;
         inter = zeros(25,1);
        for ii = 1:5
            for jj = 1:5
         
                F11 = x1(1)*dphat1dz1(gauss5(ii),gauss5(jj)) + x1(2)*dphat2dz1(gauss5(ii),gauss5(jj)) + x1(3)*dphat3dz1(gauss5(ii),gauss5(jj)) + x1(4)*dphat4dz1(gauss5(ii),gauss5(jj));
                F12 = x1(1)*dphat1dz2(gauss5(ii),gauss5(jj)) + x1(2)*dphat2dz2(gauss5(ii),gauss5(jj)) + x1(3)*dphat3dz2(gauss5(ii),gauss5(jj)) + x1(4)*dphat4dz2(gauss5(ii),gauss5(jj));
                F21 = x2(1)*dphat1dz1(gauss5(ii),gauss5(jj)) + x2(2)*dphat2dz1(gauss5(ii),gauss5(jj)) + x2(3)*dphat3dz1(gauss5(ii),gauss5(jj)) + x2(4)*dphat4dz1(gauss5(ii),gauss5(jj));
                F22 = x2(1)*dphat1dz2(gauss5(ii),gauss5(jj)) + x2(2)*dphat2dz2(gauss5(ii),gauss5(jj)) + x2(3)*dphat3dz2(gauss5(ii),gauss5(jj)) + x2(4)*dphat4dz2(gauss5(ii),gauss5(jj));
                
                F = [F11 F12;
                     F21 F22];
                
                 dphatdz1vec = [dphat1dz1(gauss5(ii),gauss5(jj)) dphat2dz1(gauss5(ii),gauss5(jj)) dphat3dz1(gauss5(ii),gauss5(jj)) dphat4dz1(gauss5(ii),gauss5(jj))];
                 dphatdz2vec = [dphat1dz2(gauss5(ii),gauss5(jj)) dphat2dz2(gauss5(ii),gauss5(jj)) dphat3dz2(gauss5(ii),gauss5(jj)) dphat4dz2(gauss5(ii),gauss5(jj))];
                 
                 phatvec = [phat1(gauss5(ii),gauss5(jj)) phat2(gauss5(ii),gauss5(jj)) phat3(gauss5(ii),gauss5(jj)) phat4(gauss5(ii),gauss5(jj))];
                
                 
                 inter1 = (1/deltat)*phatvec(A)*rho;
                 inter2 = cp*phatvec(B);
                 
                 
                 inter(count) = w5(ii)*w5(jj)*inter1*inter2*(det(F));
                 
         fM1(A,B) = sum(inter);
           
           count = count+1;
           
            end
        end
     end
 end

 Mstore1(e) = {fM1};
 
M(conn(e,:),conn(e,:)) = M(conn(e,:),conn(e,:)) + fM1;
end


%Boundary Conditions for M
for e=1:Nr:Ne

% 
 M(conn(e,1),:) = 0;
 
 M(conn(e,4),:) = 0;
 
 M(conn(e,1),conn(e,1)) = 1;
 
 M(conn(e,4),conn(e,4)) = 1;
 
 
 end






%Lumped Mass Matrix
for dia = 1:1:length(M)
    
   M2(dia,dia) = sum(M(dia,:)); 
     
    
end

%Lumped Mass Inverse

M22 = inv(M2);



%% Solve Initial FEM Matrix 

a = Tbar*ones(NP,1);

%a(1:Ne) = Tbar;

%Forward Euler Timestepping

for t=1:1:steps
    
   a(:,t+1) = M22*((M2-K)*a(:,t) + R);
   
   %a(1:Nr:Ne,t+1) = Tbar;
    
    
    
end



%Plot the solution
figure(1)
patch('Vertices', [nodecart a(:,1)],'Faces', conn, 'FaceVertexCData', a(:,1), 'FaceCol', 'interp')
colorbar

figure(2)
patch('Vertices', [nodecart a(:,101)],'Faces', conn, 'FaceVertexCData', a(:,101), 'FaceCol', 'interp')
colorbar


figure(3)
patch('Vertices', [nodecart a(:,201)],'Faces', conn, 'FaceVertexCData', a(:,201), 'FaceCol', 'interp')
colorbar


figure(4)
patch('Vertices', [nodecart a(:,301)],'Faces', conn, 'FaceVertexCData', a(:,301), 'FaceCol', 'interp')
colorbar

figure(5)
patch('Vertices', [nodecart a(:,401)],'Faces', conn, 'FaceVertexCData', a(:,401), 'FaceCol', 'interp')
colorbar

figure(6)
patch('Vertices', [nodecart a(:,501)],'Faces', conn, 'FaceVertexCData', a(:,501), 'FaceCol', 'interp')
colorbar

figure(7)
patch('Vertices', [nodecart a(:,end)],'Faces', conn, 'FaceVertexCData', a(:,end), 'FaceCol', 'interp')
colorbar



%Cross-section
rrr = nodepolar(1:Nr+1,1);
Tfem1 = a(1:Nr+1,1);
Tfem2 = a(1:Nr+1,101);
Tfem3 = a(1:Nr+1,201);
Tfem4 = a(1:Nr+1,301);
Tfem5 = a(1:Nr+1,401);
Tfem6 = a(1:Nr+1,501);
Tfem7 = a(1:Nr+1,601);
Tfem8 = a(1:Nr+1,1001);
Tfem9 = a(1:Nr+1,5001);
Tfem10 = a(1:Nr+1,end);
% 

figure(8)

hold all
plot(rrr,Tfem1)
plot(rrr,Tfem2)
plot(rrr,Tfem3)
plot(rrr,Tfem4)
plot(rrr,Tfem5)
plot(rrr,Tfem6)
plot(rrr,Tfem7)
plot(rrr,Tfem8)
plot(rrr,Tfem9)
plot(rrr,Tfem10,'k')

xlabel('Radius')
ylabel('Temperature')
legend('t = 0','t = 1','t = 2','t = 3','t = 4','t = 5','t = 6','t = 10','t = 50','t = 100')
