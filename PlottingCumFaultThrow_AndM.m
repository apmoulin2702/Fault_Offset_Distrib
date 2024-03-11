%% The purpose of this script is to plot the throw of individual faults
%% and their cumulative sum along the axis of the Red Sea


% draw lines along base and top of individual faults in QGIS and ascribe ID
% to each fault
% for each ID, compute the vertical throw (T) and distance from rotation
% pole (D)
% obtain a 2 columns matrix containing ID, D and T



load 'BASE_flankE.txt'
FB=BASE_flankE; % 4 columns (fault ID, LONG, LAT, depth)

load 'TOP_flankE.txt'
FT=TOP_flankE; % 4 columns (fault ID, LONG, LAT, depth)

load 'BaseCode_flankE.txt'
BC=BaseCode_flankE;

load 'TopCode_flankE.txt'
TC=TopCode_flankE;

for i=1:length(FB)
    [M,I]=min(abs(BC(:,2)-(FB(i,1)+1)));
    FB(i,5)=BC(I,1); % store corrected fault ID
    X1(i)=M;
end

for i=1:length(FT)
    [M,I]=min(abs(TC(:,2)-(FT(i,1)+1)));
    FT(i,5)=TC(I,1); % store corrected fault ID
    X2(i)=M;
end

load 'RedSea_SegmentedAxis.txt'
RSSA=RedSea_SegmentedAxis;


%% Compute longitude and latitude in pole-centered coordinate system

% Prepare rotation coordinates

LONG_P=deg2rad(24.22);    % Modify rotation parameters of Nubia-Arabia
LAT_P=deg2rad(31.61);     % Modify rotation parameters of Nubia-Arabia
omega_P=deg2rad(0.387);

x_MOD=cos(LAT_P).*cos(LONG_P);        % Conversion of rotation pole to cartesian
y_MOD=cos(LAT_P).*sin(LONG_P);
z_MOD=sin(LAT_P);

x_MOD_2=y_MOD/sqrt(y_MOD^2+x_MOD^2);        % Calculate cartesian coordinates of the pole about which the data have to be rotated in order for them to be "centered on the Arabia/Nubia rotation pole"
y_MOD_2=-x_MOD/sqrt(y_MOD^2+x_MOD^2);
z_MOD_2=0;


RR=deg2rad(-90)+LAT_P;         % Define rotation matrix RM_MOD
c=cos(RR);
s=sin(RR);

RM_MOD=[x_MOD_2^2*(1-c)+c x_MOD_2*y_MOD_2*(1-c)-z_MOD_2*s x_MOD_2*z_MOD_2*(1-c)+y_MOD_2*s;x_MOD_2*y_MOD_2*(1-c)+z_MOD_2*s y_MOD_2^2*(1-c)+c y_MOD_2*z_MOD_2*(1-c)-x_MOD_2*s;x_MOD_2*z_MOD_2*(1-c)-y_MOD_2*s y_MOD_2*z_MOD_2*(1-c)+x_MOD_2*s z_MOD_2^2*(1-c)+c];

% Rotate coordinates

FT(:,6)=deg2rad(FT(:,2));        % Conversion LONG in radians
FT(:,7)=deg2rad(FT(:,3));        % Conversion LAT in radians
FB(:,6)=deg2rad(FB(:,2));        % Conversion LONG in radians
FB(:,7)=deg2rad(FB(:,3));        % Conversion LAT in radians
RSSA(:,4)=deg2rad(RSSA(:,1));        % Conversion LONG in radians
RSSA(:,5)=deg2rad(RSSA(:,2));        % Conversion LAT in radians

R=6400000;
x_cart=R*cos(FT(:,7)).*cos(FT(:,6));        % Conversion of LONG and LAT to cartesian
y_cart=R*cos(FT(:,7)).*sin(FT(:,6));
z_cart=R*sin(FT(:,7));
V=horzcat(x_cart,y_cart,z_cart);        % Store cartesian coordinates in vector V

V_MOD=V*RM_MOD;        % Rotate cartesian coordinates
[LONG_MOD,LAT_MOD,R_MOD]=cart2sph(V_MOD(:,1),V_MOD(:,2),V_MOD(:,3));        % Convert coordinates to geographic

FT(:,8)=rad2deg(LONG_MOD);         % Store LONG relative to rotation pole
FT(:,9)=rad2deg(LAT_MOD);         % Store LAT relative to rotation pole (opposite to distance from rotation pole)

clear LONG_MOD
clear LAT_MOD
clear R_MOD
clear V_MOD
clear x_cart
clear y_cart
clear z_cart
clear V

x_cart=R*cos(FB(:,7)).*cos(FB(:,6));        % Conversion of LONG and LAT to cartesian
y_cart=R*cos(FB(:,7)).*sin(FB(:,6));
z_cart=R*sin(FB(:,7));
V=horzcat(x_cart,y_cart,z_cart);        % Store cartesian coordinates in vector V

V_MOD=V*RM_MOD;        % Rotate cartesian coordinates
[LONG_MOD,LAT_MOD,R_MOD]=cart2sph(V_MOD(:,1),V_MOD(:,2),V_MOD(:,3));        % Convert coordinates to geographic

FB(:,8)=rad2deg(LONG_MOD);         % Store LONG relative to rotation pole
FB(:,9)=rad2deg(LAT_MOD);         % Store LAT relative to rotation pole (opposite to distance from rotation pole)

clear LONG_MOD
clear LAT_MOD
clear R_MOD
clear V_MOD
clear x_cart
clear y_cart
clear z_cart
clear V

x_cart=R*cos(RSSA(:,5)).*cos(RSSA(:,4));        % Conversion of LONG and LAT to cartesian
y_cart=R*cos(RSSA(:,5)).*sin(RSSA(:,4));
z_cart=R*sin(RSSA(:,5));
V=horzcat(x_cart,y_cart,z_cart);        % Store cartesian coordinates in vector V

V_MOD=V*RM_MOD;        % Rotate cartesian coordinates
[LONG_MOD,LAT_MOD,R_MOD]=cart2sph(V_MOD(:,1),V_MOD(:,2),V_MOD(:,3));        % Convert coordinates to geographic

RSSA(:,7)=rad2deg(LONG_MOD);         % Store LONG relative to rotation pole
RSSA(:,8)=rad2deg(LAT_MOD);         % Store LAT relative to rotation pole (opposite to distance from rotation pole)


%% load IZ data

load ITZ.txt
ITZ=ITZ;

% Rotate coordinates

LONG_ITZ=deg2rad(ITZ(:,1));        % Conversion LONG in radians
LAT_ITZ=deg2rad(ITZ(:,2));        % Conversion LAT in radians

R=6400000;
x_cart_itz=R*cos(LAT_ITZ).*cos(LONG_ITZ);        % Conversion of LONG and LAT to cartesian
y_cart_itz=R*cos(LAT_ITZ).*sin(LONG_ITZ);
z_cart_itz=R*sin(LAT_ITZ);
V_ITZ=horzcat(x_cart_itz,y_cart_itz,z_cart_itz);        % Store cartesian coordinates in vector V_ITZ


V_ITZ_MOD=V_ITZ*RM_MOD;        % Rotate cartesian coordinates
[LONG_ITZ_MOD,LAT_ITZ_MOD,R_ITZ_MOD]=cart2sph(V_ITZ_MOD(:,1),V_ITZ_MOD(:,2),V_ITZ_MOD(:,3));        % Convert coordinates to geographic

LONG_ITZ_MOD=rad2deg(LONG_ITZ_MOD);         % Store LONG relative to rotation pole
LAT_ITZ_MOD=rad2deg(LAT_ITZ_MOD);         % Store LAT relative to rotation pole (opposite to distance from rotation pole)

ITZ(:,3)=LONG_ITZ_MOD;
ITZ(:,4)=LAT_ITZ_MOD;

DIST_ITZ=LONG_ITZ_MOD.*((2*pi*R.*cos(deg2rad(LAT_ITZ_MOD)))/360); % compute distance in km along rotation pole-centered parallels

ITZ(:,5)=DIST_ITZ;



%% Extract throw data for each fault

for i=min(FB(:,5)):max(FB(:,5))

    clear M1
    clear M2

    M1=FB;
    M1(M1(:,5)<i,:)=[]; 
    M1(M1(:,5)>i,:)=[]; % extract sub-matrix containing the data of the base of fault "i" 

    M2=FT;
    M2(M2(:,5)<i,:)=[];
    M2(M2(:,5)>i,:)=[]; % extract sub-matrix containing the data of the top of fault "i" 

VT{i}=zeros(size(M2,1),4);

for j=1:size(M2,1)
    [~,ind]=min(abs(M1(:,9)-M2(j,9))); % find point in "base" closest (in terms of rotation-pole centered latitude) to targeted point in "top"
    %% put a threshold for maximum angular distance between top and base
    VT{i}(j,1)=M2(j,9); % store rotation-pole centered latitude
    VT{i}(j,2)=M2(j,8); % store rotation-pole centered longitude
    VT{i}(j,3)=i; % store fault ID
    VT{i}(j,4)=M2(j,4)-M1(ind,4); % compute vertical separation between "top" and "base" at the targeted point
end

end


%% plot throw of individual faults as a function of distance from rotation pole

figure;

for i=min(FB(:,5)):max(FB(:,5))
    plot(VT{i}(:,1),VT{i}(:,4),'k')
    hold on
end


%% concatenate data from all faults

VTcat=VT{min(FB(:,5))};

for i=min(FB(:,5)):max(FB(:,5))-1
    VTcat=vertcat(VTcat,VT{i+1});
end


%% calculate cumulative throw against latitude by analyzing throw distribution within a moving window

W=0.005; % modify sample spacing if needed (in °)
Width=0.005; % modify with of averaging window if needed (in °)
LAT1=min(VTcat(:,1))-W; 
D_LAT=max(VTcat(:,1))-min(VTcat(:,1)); % compute range of latitude

Mat=zeros(0,5); % define matrix where results will be stored
Temp2=zeros(0,5);

for j=1:floor(D_LAT/W)

clear Temp
LAT1=LAT1+W; % define lower bound of LAT
LAT2=LAT1+Width; % define upper bound of LAT
Temp=VTcat(VTcat(:,1)>LAT1 & VTcat(:,1)<LAT2,:); % extract data corresponding to targeted interval
UQ=unique(Temp(:,3)); % store the list of faults ID found in targeted interval
clear Temp2
clear M_UQ

for i=1:length(UQ)
    
    M_UQ{i}=Temp(Temp(:,3)==UQ(i),:); % generate one sub-matrix per fault
    
    Temp2(i,1)=(LAT1+LAT2)/2; % compute average distance from pole
    Temp2(i,2)=mean(M_UQ{i}(:,2)); % compute average distance from ridge
    Temp2(i,3)=mean(M_UQ{i}(:,4)); % compute average vertical throw

end

if exist('Temp2')

Temp2=sortrows(Temp2,3); % sort faults as a function of their distance to the ridge

for i=1:length(UQ)
    Temp2(i,4)=i; % define order of fault (increasing from 1 away from the ridge)
    Temp2(i,5)=sum(Temp2(1:i,3)); % calculate cumulative throw
end

Mat=vertcat(Mat,Temp2);

end

end


%% complete matrix with missing values

UQ2=unique(Mat(:,4)); % determine maximum number of fault orders in parallel and store the list of fault orders 
x=min(VTcat(:,1))+Width/2:W:min(VTcat(:,1))+W*floor(D_LAT/W)'; % define the vector containing sampling latitudes

for i=1:length(UQ2)

Temp4=Mat(Mat(:,4)==UQ2(i),:); % extract sub-matrix of the targeted fault order
[val,pos]=setdiff(round(x,4),round(Temp4(:,1),4)); % extract values and position of the latitudes that are not available for this fault order

Mat_Supp{i}(:,1)=val; % store missing latitude values in first column of Mat_Supp
Mat_Supp{i}(:,2)=NaN(length(val),1);
Mat_Supp{i}(:,3)=NaN(length(val),1);
Mat_Supp{i}(:,4)=repmat(UQ2(i),length(val),1);
Mat_Supp{i}(:,5)=zeros(length(val),1); 

Mat=vertcat(Mat,Mat_Supp{i}); % add missing value to original matrix

end

Mat(:,6)=Mat(:,2).*cos(deg2rad(90-Mat(:,1)))*111.325; % compute and store relative length along parallel (in km)


%% plot cumulative throw profiles

UQ3=unique(Mat(:,4)); % store the list of fault order in final data-set

figure;
hold on

r=ones(length(UQ3)+1,1); % define new color scale
g=linspace(0,1,length(UQ3)+1); % define new color scale (continuing)
b=g; % define new color scale (continuing)
redmap=[r g' b']; % define new color scale (continuing)

for i=1:length(UQ3)
    clear Temp3
    Temp3=Mat(Mat(:,4)==UQ3(length(UQ3)+1-i),:); % extract sub-matrix with unique fault order (in descending order)
    area(Temp3(:,1),Temp3(:,5),'FaceColor',redmap(i,:)); % color are under the curve with above-defined colorscale varying with distance from the ridge
end


%% plot all throw together with averaged max throw

% extract averaged max throw along strike

XX=min(VTcat(:,1)):0.02:min(VTcat(:,1))+0.02*(round((max(VTcat(:,1))-min(VTcat(:,1)))/0.02));

for i=1:length(XX)
    
    clear Temp8
    Temp8=VTcat(VTcat(:,1)>XX(i)-0.05 & VTcat(:,1)<XX(i)+0.05,:);
    
    if size(Temp8,1)>0
    max_T(i)=max(Temp8(:,4));
    else
        max_T(i)=NaN;
    end

end


% plot color-coded area under curved together with individual throws

max_T2=max_T*5;
%% Note that max throw has been multiplied by 5 to get a scale equivalent to individual throws data
max_T2(isnan(max_T2))=0;

n_vec=min(XX):0.02:max(XX); % define resolution of grid
x_vec = zeros(1,length(XX)); % define bases (intersection with y axis) of y values
y_vec = max_T2;
N=[n_vec,fliplr(n_vec)];
X=[x_vec,flip(y_vec)];

y=[zeros(1,length(max_T2))';max_T2'];
resolution=[2000,2000];
px=linspace(min(n_vec),max(n_vec),resolution(1));
py_=linspace(min(y),max(y),resolution(2));
[px,py_]=meshgrid(px,py_);

in=inpolygon(px,py_,N,X);

pz=py_;
pz(~in)=NaN;


%% Quantify magmatic contribution
   
% compute normalized cumulative throw

W=0.005; % modify sample spacing if needed
Width=0.005; % modify width of averaging window if needed
LAT1=min(VTcat(:,1))-W; 
D_LAT=max(VTcat(:,1))-min(VTcat(:,1)); % compute range of latitude

for j=1:floor(D_LAT/W) % for each targeted interval

clear Temp9
clear Temp10
clear maxID
clear hDw
clear maxThrow
clear AX_DIST
clear AX_DIST_2
LAT1=LAT1+W; % define lower bound of LAT
LAT2=LAT1+Width; % define upper bound of LAT
Temp9=Mat(Mat(:,1)>LAT1 & Mat(:,1)<LAT2,:); % extract data corresponding to targeted interval
Temp9(Temp9(:,5)==0,:)=[]; % keep data with non-null cumulative throw
if length(Temp9)>0
maxID=max(Temp9(:,4)); % extract maximum ID found in targeted interval
Temp10=Temp9(Temp9(:,4)==maxID,:); % extract row corresponding to maxID in targeted interval
maxThrow=max(Temp10(:,5)); % extract max cumluative throw in this data-set
Temp10=Temp10(Temp10(:,5)==maxThrow,:); % keep row with max cumulative throw
    if size(Temp10,1)>1
        Temp10(2,:)=[];
    end
    [dist,IND]=min(abs(RSSA(:,8)-Temp10(:,1))); % extract point along the RS axis which has a rotation-pole-centered latitude which is closest to that of the point of interest
    AX_DIST=abs(RSSA(IND,7)-Temp10(:,2)); % compute distance from ridge of the point of interest (in ° around a small circle)
    AX_DIST_2=(deg2rad(AX_DIST)/(2*pi))*2*pi*R*cos(deg2rad(Temp10(:,1))); % compute distance from ridge of the point of interest (in m around a small circle)
norm_throw_e(j,1)=Temp10(:,1); % store latitude from rotation pole
norm_throw_e(j,2)=Temp10(:,5)/AX_DIST_2; % compute normalized throw as the cumulative throw recorded by fault of largest ID divided by the opening recorded by this fault
norm_throw_e(j,3)=1-norm_throw_e(j,2)/cos(deg2rad(60)); % compute magmatic contribution (M) to opening assuming fault dips of 60°
norm_throw_e(j,4)=AX_DIST_2; % store distance from ridge
end

end


% remove normalized cumulative throw within ITZ

for i=1:length(ITZ)/2
    minB=ITZ(1+2*(i-1),4); % extract lat of southern ITZ limit
    maxB=ITZ((1+2*(i-1))+1,4); % extract lat of northern ITZ limit
    Temp11=norm_throw_e(norm_throw_e(:,1)<minB,:); % keep data south of minB
    Temp12=norm_throw_e(norm_throw_e(:,1)>maxB,:); % keep data north of maxB
    norm_throw_e=vertcat(Temp11,Temp12); % store data outside of this ITZ
end

norm_throw_e(norm_throw_e(:,1)==0,:)=[];


% remove sampled points that are too close to the ridge

r=1;
for j=1:length(norm_throw_e)
    if norm_throw_e(j,4)>6000 % if sampled point is farther from the ridge than a threshold (in m)
        M_e(r,1)=norm_throw_e(j,1); % store latitude
        M_e(r,2)=norm_throw_e(j,3); % store magmatic contribution M
        r=r+1;
    end
end


% averaged magmatic contribution (M)

window=0.05; % modify averaging window if needed
LAT1=min(VTcat(:,1))-window/2;
LAT2=LAT1+window;
D_LAT=max(VTcat(:,1))-min(VTcat(:,1)); % compute range of latitude

for j=1:2*floor(D_LAT/window)
    clear Temp13
    LAT1=LAT1+window/2;
    LAT2=LAT1+window;
    Temp13=M_e(M_e(:,1)>LAT1 & M_e(:,1)<LAT2,:); % extract data corresponding to targeted interval
    avM_e(j,1)=(LAT1+LAT2)/2;
    avM_e(j,2)=mean(Temp13(:,2)); % extract average of M within the targeted interval
end



%% Store data 

XX1=px;
YY1=py_;
ZZ1=pz;
L1=max_T2;
xXX=XX;
Mat1=Mat;

clearvars -except XX1 YY1 ZZ1 L1 xXX Mat1 M_e avM_e RSSA ITZ



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% SAME WITH WESTERN FLANK %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load 'BASE_flankW.txt'
FB=BASE_flankW; % 4 columns (fault ID, LONG, LAT, depth)

load 'TOP_flankW.txt'
FT=TOP_flankW; % 4 columns (fault ID, LONG, LAT, depth)


%% Compute longitude and latitude in pole-centered coordinate system

% Prepare rotation coordinates

LONG_P=deg2rad(24.22);    % Modify rotation parameters of Nubia-Arabia
LAT_P=deg2rad(31.61);     % Modify rotation parameters of Nubia-Arabia
omega_P=deg2rad(0.387);

x_MOD=cos(LAT_P).*cos(LONG_P);        % Conversion of rotation pole to cartesian
y_MOD=cos(LAT_P).*sin(LONG_P);
z_MOD=sin(LAT_P);

x_MOD_2=y_MOD/sqrt(y_MOD^2+x_MOD^2);        % Calculate cartesian coordinates of the pole about which the data have to be rotated in order for them to be "centered on the Arabia/Nubia rotation pole"
y_MOD_2=-x_MOD/sqrt(y_MOD^2+x_MOD^2);
z_MOD_2=0;


RR=deg2rad(-90)+LAT_P;         % Define rotation matrix RM_MOD
c=cos(RR);
s=sin(RR);

RM_MOD=[x_MOD_2^2*(1-c)+c x_MOD_2*y_MOD_2*(1-c)-z_MOD_2*s x_MOD_2*z_MOD_2*(1-c)+y_MOD_2*s;x_MOD_2*y_MOD_2*(1-c)+z_MOD_2*s y_MOD_2^2*(1-c)+c y_MOD_2*z_MOD_2*(1-c)-x_MOD_2*s;x_MOD_2*z_MOD_2*(1-c)-y_MOD_2*s y_MOD_2*z_MOD_2*(1-c)+x_MOD_2*s z_MOD_2^2*(1-c)+c];

% Rotate coordinates

FT(:,6)=deg2rad(FT(:,2));        % Conversion LONG in radians
FT(:,7)=deg2rad(FT(:,3));        % Conversion LAT in radians
FB(:,6)=deg2rad(FB(:,2));        % Conversion LONG in radians
FB(:,7)=deg2rad(FB(:,3));        % Conversion LAT in radians

R=6400000;
x_cart=R*cos(FT(:,7)).*cos(FT(:,6));        % Conversion of LONG and LAT to cartesian
y_cart=R*cos(FT(:,7)).*sin(FT(:,6));
z_cart=R*sin(FT(:,7));
V=horzcat(x_cart,y_cart,z_cart);        % Store cartesian coordinates in vector V

V_MOD=V*RM_MOD;        % Rotate cartesian coordinates
[LONG_MOD,LAT_MOD,R_MOD]=cart2sph(V_MOD(:,1),V_MOD(:,2),V_MOD(:,3));        % Convert coordinates to geographic

FT(:,8)=rad2deg(LONG_MOD);         % Store LONG relative to rotation pole
FT(:,9)=rad2deg(LAT_MOD);         % Store LAT relative to rotation pole (opposite to distance from rotation pole)

clear LONG_MOD
clear LAT_MOD
clear R_MOD
clear V_MOD
clear x_cart
clear y_cart
clear z_cart
clear V

x_cart=R*cos(FB(:,7)).*cos(FB(:,6));        % Conversion of LONG and LAT to cartesian
y_cart=R*cos(FB(:,7)).*sin(FB(:,6));
z_cart=R*sin(FB(:,7));
V=horzcat(x_cart,y_cart,z_cart);        % Store cartesian coordinates in vector V

V_MOD=V*RM_MOD;        % Rotate cartesian coordinates
[LONG_MOD,LAT_MOD,R_MOD]=cart2sph(V_MOD(:,1),V_MOD(:,2),V_MOD(:,3));        % Convert coordinates to geographic

FB(:,8)=rad2deg(LONG_MOD);         % Store LONG relative to rotation pole
FB(:,9)=rad2deg(LAT_MOD);         % Store LAT relative to rotation pole (opposite to distance from rotation pole)


%% Extract throw data for each fault

for i=min(FB(:,1)):max(FB(:,1))

    clear M1
    clear M2

    M1=FB;
    M1(M1(:,1)<i,:)=[]; 
    M1(M1(:,1)>i,:)=[]; % extract sub-matrix containing the data of the base of fault "i" 

    M2=FT;
    M2(M2(:,1)<i,:)=[];
    M2(M2(:,1)>i,:)=[]; % extract sub-matrix containing the data of the top of fault "i" 

VT{i}=zeros(size(M2,1),4);

for j=1:size(M2,1)
    [~,ind]=min(abs(M1(:,9)-M2(j,9))); % find point in "base" closest (in terms of rotation-pole centered latitude) to targeted point in "top"
    %% put a threshold for maximum angular distance between top and base
    VT{i}(j,1)=M2(j,9); % store rotation-pole centered latitude
    VT{i}(j,2)=M2(j,8); % store rotation-pole centered longitude
    VT{i}(j,3)=i; % store fault ID
    VT{i}(j,4)=M2(j,4)-M1(ind,4); % compute vertical separation between "top" and "base" at the targeted point
end

end


%% plot throw of individual faults as a function of distance from rotation pole

figure;

for i=min(FB(:,1)):max(FB(:,1))
    plot(VT{i}(:,1),VT{i}(:,4),'k')
    hold on
end


%% concatenate data from all faults

VTcat=VT{min(FB(:,1))};

for i=min(FB(:,1)):max(FB(:,1))-1
    VTcat=vertcat(VTcat,VT{i+1});
end


%% calculate cumulative throw against latitude by analyzing throw distribution within a moving window

W=0.005; % modify sample spacing if needed
Width=0.02; % modify width of averaging window if needed
LAT1=min(VTcat(:,1))-W; 
D_LAT=max(VTcat(:,1))-min(VTcat(:,1)); % compute range of latitude

Mat=zeros(0,5); % define matrix where results will be stored
Temp2=zeros(0,5);

for j=1:floor(D_LAT/W)

clear Temp
LAT1=LAT1+W; % define lower bound of LAT
LAT2=LAT1+Width; % define upper bound of LAT
Temp=VTcat(VTcat(:,1)>LAT1 & VTcat(:,1)<LAT2,:); % extract data corresponding to targeted interval
UQ=unique(Temp(:,3)); % store the list of faults ID found in targeted interval
clear Temp2
clear M_UQ


for i=1:length(UQ)
    
    M_UQ{i}=Temp(Temp(:,3)==UQ(i),:); % generate one sub-matrix per fault
    
    Temp2(i,1)=(LAT1+LAT2)/2; % compute average distance from pole
    Temp2(i,2)=mean(M_UQ{i}(:,2)); % compute average distance from ridge
    Temp2(i,3)=mean(M_UQ{i}(:,4)); % compute average vertical throw

end


if exist('Temp2')

Temp2=sortrows(Temp2,3); % sort faults as a function of their distance to the ridge

for i=1:length(UQ)
    Temp2(i,4)=i; % define order of fault (increasing from 1 away from the ridge)
    Temp2(i,5)=sum(Temp2(1:i,3)); % calculate cumulative throw
end

Mat=vertcat(Mat,Temp2);

end

end


%% complete matrix with missing values

UQ2=unique(Mat(:,4)); % determine maximum number of fault orders in parallel and store the list of fault orders 
x=min(VTcat(:,1))+Width/2:W:min(VTcat(:,1))+W*floor(D_LAT/W)'; % define the vector containing sampling latitudes

for i=1:length(UQ2)

Temp4=Mat(Mat(:,4)==UQ2(i),:); % extract sub-matrix of the targeted fault order
[val,pos]=setdiff(round(x,4),round(Temp4(:,1),4)); % extract values and position of the latitudes that are not available for this fault order

Mat_Supp{i}(:,1)=val; % store missing latitude values in first column of Mat_Supp
Mat_Supp{i}(:,2)=NaN(length(val),1);
Mat_Supp{i}(:,3)=NaN(length(val),1);
Mat_Supp{i}(:,4)=repmat(UQ2(i),length(val),1);
Mat_Supp{i}(:,5)=zeros(length(val),1); 

Mat=vertcat(Mat,Mat_Supp{i}); % add missing value to original matrix

end

Mat(:,6)=Mat(:,2).*cos(deg2rad(90-Mat(:,1)))*111.325; % compute and store relative length along parallel (in km)


%% plot cumulative throw profiles

UQ3=unique(Mat(:,4)); % store the list of fault order in final data-set

figure;
hold on

r=ones(length(UQ3)+1,1); % define new color scale
g=linspace(0,1,length(UQ3)+1); % define new color scale (continuing)
b=g; % define new color scale (continuing)
redmap=[r g' b']; % define new color scale (continuing)

for i=1:length(UQ3)
    clear Temp3
    Temp3=Mat(Mat(:,4)==UQ3(length(UQ3)+1-i),:); % extract sub-matrix with unique fault order (in descending order)
    area(Temp3(:,1),Temp3(:,5),'FaceColor',redmap(i,:)); % color are under the curve with above-defined colorscale varying with distance from the ridge
end


%% Quantify magmatic contribution

% compute normalized cumulative throw

W=0.005; % modify sample spacing if needed
Width=0.005; % modify width of averaging window if needed
LAT1=min(VTcat(:,1))-W; 
D_LAT=max(VTcat(:,1))-min(VTcat(:,1)); % compute range of latitude

for j=1:floor(D_LAT/W) % for each targeted interval

clear Temp9
clear Temp10
clear maxID
clear hDw
clear maxThrow
LAT1=LAT1+W; % define lower bound of LAT
LAT2=LAT1+Width; % define upper bound of LAT
Temp9=Mat(Mat(:,1)>LAT1 & Mat(:,1)<LAT2,:); % extract data corresponding to targeted interval
Temp9(Temp9(:,5)==0,:)=[]; % keep data with non-null cumulative throw
if length(Temp9)>0
maxID=max(Temp9(:,4)); % extract maximum ID found in targeted interval
Temp10=Temp9(Temp9(:,4)==maxID,:); % extract row corresponding to maxID in targeted interval
maxThrow=max(Temp10(:,5)); % extract max cumluative throw in this data-set
Temp10=Temp10(Temp10(:,5)==maxThrow,:); % keep row with max cumulative throw
    if size(Temp10,1)>1
        Temp10(2,:)=[];
    end
    [dist,IND]=min(abs(RSSA(:,8)-Temp10(:,1))); % extract point along the RS axis which has a rotation-pole-centered latitude which is closest to that of the point of interest
    AX_DIST=abs(RSSA(IND,7)-Temp10(:,2)); % compute distance from ridge of the point of interest (in ° around a small circle)
    AX_DIST_2=(deg2rad(AX_DIST)/(2*pi))*2*pi*R*cos(deg2rad(Temp10(:,1))); % compute distance from ridge of the point of interest (in m around a small circle)
norm_throw_w(j,1)=Temp10(:,1); % store latitude from rotation pole
norm_throw_w(j,2)=Temp10(:,5)/AX_DIST_2; % compute normalized throw as the cumulative throw recorded by fault of largest ID divided by the opening recorded by this fault
norm_throw_w(j,3)=1-norm_throw_w(j,2)/cos(deg2rad(60)); % compute magmatic contribution (M) to opening assuming fault dips of 60°
norm_throw_w(j,4)=AX_DIST_2; % store distance from ridge
end

end


% remove normalized cumulative throw within ITZ

for i=1:length(ITZ)/2
    minB=ITZ(1+2*(i-1),4); % extract lat of southern ITZ limit
    maxB=ITZ((1+2*(i-1))+1,4); % extract lat of northern ITZ limit
    Temp11=norm_throw_w(norm_throw_w(:,1)<minB,:); % keep data south of minB
    Temp12=norm_throw_w(norm_throw_w(:,1)>maxB,:); % keep data north of maxB
    norm_throw_w=vertcat(Temp11,Temp12); % store data outside of this ITZ
end

norm_throw_w(norm_throw_w(:,1)==0,:)=[];


% remove sampled points that are too close to the ridge

r=1;
for j=1:length(norm_throw_w)
    if norm_throw_w(j,4)>6000 % if sampled point is farther from the ridge than a threshold (in m)
        M_w(r,1)=norm_throw_w(j,1); % store latitude
        M_w(r,2)=norm_throw_w(j,3); % store magmatic contribution M
        r=r+1;
    end
end

% averaged magmatic contribution (M)

window=0.05; % modify averaging window if needed
LAT1=min(VTcat(:,1))-window/2;
LAT2=LAT1+window;
D_LAT=max(VTcat(:,1))-min(VTcat(:,1)); % compute range of latitude

for j=1:2*floor(D_LAT/window)
    clear Temp13
    LAT1=LAT1+window/2;
    LAT2=LAT1+window;
    Temp13=M_w(M_w(:,1)>LAT1 & M_w(:,1)<LAT2,:); % extract data corresponding to targeted interval
    avM_w(j,1)=(LAT1+LAT2)/2;
    avM_w(j,2)=mean(Temp13(:,2)); % extract average of M within the targeted interval
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% plot all throw together with averaged max throw

% extract averaged max throw along strike

XX=min(VTcat(:,1)):0.02:min(VTcat(:,1))+0.02*(round((max(VTcat(:,1))-min(VTcat(:,1)))/0.02));

for i=1:length(XX)
    
    clear Temp8
    Temp8=VTcat(VTcat(:,1)>XX(i)-0.05 & VTcat(:,1)<XX(i)+0.05,:); % extract dataset within latitude interval
    
    if size(Temp8,1)>0
    max_T(i)=max(Temp8(:,4)); % extract maximum throw
    else
        max_T(i)=NaN;
    end

end


% plot color-coded area under curved together with individual throws

max_T2=max_T*5; % multiply maximum throw (for scaling purpose)
%% Note that max throw has been multiplied by 5 to get a scale equivalent to individual throws data
max_T2(isnan(max_T2))=0;

n_vec=min(XX):0.02:max(XX); % define resolution of grid
x_vec = zeros(1,length(XX)); % define bases (intersection with y axis) of y values
y_vec = max_T2; 
N=[n_vec,fliplr(n_vec)]; 
X=[x_vec,flip(y_vec)];

y=[zeros(1,length(max_T2))';max_T2'];
resolution=[2000,2000];
px=linspace(min(n_vec),max(n_vec),resolution(1));
py_=linspace(min(y),max(y),resolution(2));
[px,py_]=meshgrid(px,py_);

in=inpolygon(px,py_,N,X);

pz=py_;
pz(~in)=NaN;

figure;
surf(px,py_,pz,'edgecolor','none');
hold on
plot3(XX,max_T2,max_T2,'LineWidth',2,'Color','k')
view(2)


% add cumulative throw

UQ3=unique(Mat(:,4)); % store the list of fault order in final data-set

for i=1:length(UQ3)
    clear Temp3
    Temp3=Mat(Mat(:,4)==UQ3(length(UQ3)+1-i),:); % extract sub-matrix with unique fault order (in descending order)
    scatter3(Temp3(:,1),Temp3(:,5),Temp3(:,5),3,'k','filled'); % scatter cumulative throws
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% OVERLAP DATA OF THE EASTERN FLANK %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot

surf(XX1,-YY1,ZZ1,'edgecolor','none');
hold on
plot3(xXX,-L1,L1,'LineWidth',2,'Color','k')
view(2)

clear UQ3
UQ3=unique(Mat1(:,4)); % store the list of fault order in final data-set

for i=1:length(UQ3)
    clear Temp3
    Temp3=Mat1(Mat1(:,4)==UQ3(length(UQ3)+1-i),:); % extract sub-matrix with unique fault order (in descending order)
    scatter3(Temp3(:,1),-Temp3(:,5),Temp3(:,5),3,'k','filled'); % scatter cumulative throws
end



xlabel('Latitude (°) in rotation-pole-centered coordinate system','FontSize',16,'FontWeight','bold')
ylabel({'Cumulative fault offset';'computed across a 0.02°-wide';'moving window'},'FontSize',16,'FontWeight','bold')
set(gca, 'XDir','reverse')
colorbar
col = colorbar('XTickLabel', {'0', '200', '400', '600', '800', '1000', '1200', '1400'}, 'XTick',0:1000:7000);
ylabel(col,{'Maximum Fault offset (m)';'computed across a 0.1°-wide';'moving window'},'FontSize',16,'FontWeight','bold')
colormap cool
plot([68 78.5],[0 0],'black','LineWidth',1) % plot throw=0
yticklabels({'8000','6000','4000','2000','0','2000','4000','6000','8000',})


% overlap ITZ shaded bands

bands=flip(ITZ(:,4));
bands=reshape(bands,[2 12]);
bands=bands';
xb = [bands fliplr(bands)];
yb = ([[1;1]*(-8000); [1;1]*8000]*ones(1,size(bands,1))).'; 
zb = repmat(9000,12,4);

figure;

for k = 1:size(bands,1)
    patch(xb(k,:), yb(k,:), [1 1 1]*0.25, 'FaceAlpha',0.3, 'EdgeColor',[1 1 1]*0.25)
end


for k = 1:size(bands,1)
    patch(xb(k,:), yb(k,:), zb(k,:),[1 1 1]*0.25, 'FaceAlpha',0.3, 'EdgeColor',[1 1 1]*0.25)
end

% plot magmatic contribution M

figure;
scatter(M_e(:,1),1-M_e(:,2),30,'k','filled')
hold on
scatter(M_w(:,1),-1+M_w(:,2),30,'k','filled')
set(gca, 'XDir','reverse')
axis([68 78.5 -1 1])
hline = refline(0,0);
hline.Color='k';

for k = 1:size(bands,1)
    patch(xb(k,:), yb(k,:), [1 1 1]*0.25, 'FaceAlpha',0.3, 'EdgeColor',[1 1 1]*0.25) % add ITZ bands
end

xlabel('Latitude (°) in rotation-pole-centered coordinate system','FontSize',16,'FontWeight','bold')
ylabel({'Estimated magmatic';'contribution (M)'},'FontSize',16,'FontWeight','bold')

set(gca,'YTickLabel', {'0','0.5','1','0.5','0'}, 'YTick',[-1 -0.5 0 0.5 1]);





figure;
scatter3(M_e(:,1),1-M_e(:,2),repmat(9000,length(M_e),1),30,'k','filled')
hold on
scatter3(M_w(:,1),-1+M_w(:,2),repmat(9000,length(M_w),1),30,'k','filled')





