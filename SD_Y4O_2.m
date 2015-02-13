%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Indian Institute of Technology Bombay %%
%%                   CSRE                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Authors:  *Avik Bhattacharya          %%
%%           *Shaunak De                 %%
%%           *Arnab Muhuri               %%
%% Version:  *v1.0                       %%
%% Date:     *August 2014                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Call the output saved from SD_Y4O_1.m

clear all;
clc;

%% User Input

rows=input('Number of Rows: ');
cols=input('Number of Columns: ');
wsi=input('Enter Window Size: ');
look_search_range=input('Enter Look Search Range: '); % 500 or 1000
step=input('Enter step size for L: '); % 0.5 or 1

%% File Input Loop

Final=zeros(cols,rows,10);

for aa=1:1:10
    
    [filename, pathname, filterindex] = uigetfile('*.*', 'Image selection');
    fileandpath=[pathname filename];
    
    disp('   ');
    disp(['File:', fileandpath]);
    
    fid=fopen(fileandpath,'rb');
    
    Final(:,:,aa)=fread(fid,[cols rows], 'float32');
    
end

%% Output Saving Option

ch = input('Do you want to save the ouput(Y=1/N=0): ');

if(ch == 1)
    
    [filename, pathname, filterindex] = uigetfile('*.*', 'Folder to Save');
    
end

%% Image Subset

choice=input('Do you want to subset image (YES=1/NO=0): ');

if (choice == 1)
    
    row_start=input('Row Begins At: ')
    row_end=input('Row Ends At: ')
    
    col_start=input('Column Begins At: ');
    col_end=input('Column Ends At: ');
    
    rows=(row_end)-(row_start)+1
    
    cols=(col_end)-(col_start)+1
    
    Final=Final(col_start:col_end,row_start:row_end,:);
    
end

%% Image Pre-processing

ep = 0.00001;

T11=Final(:,:,1)' + ep;
T22=Final(:,:,2)' + ep;
T33=Final(:,:,3)' + ep;

ReT12=Final(:,:,4)';
ImT12=Final(:,:,5)';

T12=complex(ReT12,ImT12) + ep;
T21=conj(T12);

ReT13=Final(:,:,6)';
ImT13=Final(:,:,7)';

T13=complex(ReT13,ImT13) + ep;
T31=conj(T13);

ReT23=Final(:,:,8)';
ImT23=Final(:,:,9)';

T23=complex(ReT23,ImT23) + ep;
T32=conj(T23);

% The unwrapped angle THEETA_HD.raw file computed from SD_Y4O_1.m

Angle=Final(:,:,10)';

%% for window processing

wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopi= rows-inci; % Stop row for window processing
stopj= cols-incj; % Stop column for window processing

%% Image Processing Loop

Max_Look=zeros(rows,cols);
hd_max_look=zeros(rows,cols);

t_start = tic;
start_time = datestr(now);
h = waitbar(0,'1','Name','Progress Monitor...');
for ii=starti:stopi
    perc = (ii/(stopi))*100;
    waitbar(ii/(stopi),h, sprintf('Start Time:  %s \n%6.2f%% of Process Completed... ', start_time, perc));
    
    for jj=startj:stopj
        
        Theta = Angle(ii,jj);
        
        Bin_T11=T11(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_T11=mean(Bin_T11(:));
        
        Bin_T22=T22(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_T22=mean(Bin_T22(:));
        
        Bin_T33=T33(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_T33=mean(Bin_T33(:));
        
        %%
        Bin_T12=T12(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_T12=mean(Bin_T12(:));
        
        Bin_T21=T21(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_T21=mean(Bin_T21(:));
        
        %%
        Bin_T13=T13(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_T13=mean(Bin_T13(:));
        
        Bin_T31=T31(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_T31=mean(Bin_T31(:));
        
        %%
        Bin_T23=T23(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_T23=mean(Bin_T23(:));
        
        Bin_T32=T32(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_T32=mean(Bin_T32(:));
        
        T_Unrot= [abs(Bin_T11) Bin_T12 Bin_T13;Bin_T21 abs(Bin_T22) Bin_T23;Bin_T31 Bin_T32 abs(Bin_T33)];
        
        U = [1 0 0;0 cosd(2*Theta) sind(2*Theta);0 -sind(2*Theta) cosd(2*Theta)];
        
        T_Rot = U*T_Unrot*inv(U);
        
        T_Rot(2,2)=abs(T_Rot(2,2));
        
        T_Rot(3,3)=abs(T_Rot(3,3));
        
        %% Variation of D w.r.t. L
        
        tt = max(size(1:step:look_search_range));
        D=zeros(1,tt);
        
        zz = 0;
        for L=1:step:look_search_range
            
            zz=zz+1;
            
            D_22 = 1 - ((2*sqrt((T_Unrot(2,2))* (T_Rot(2,2))))/((T_Unrot(2,2)) + (T_Rot(2,2))))^(L);
            D_33 = 1 - ((2*sqrt((T_Unrot(3,3))* (T_Rot(3,3))))/((T_Unrot(3,3)) + (T_Rot(3,3))))^(L);
            
            temp = D_33 - D_22;
            
            if (temp < 0)
                temp = 0;
            end
            
            D(1,zz) = temp;
            
        end
        
        %% L at Maximum D
        
        [value,Index]=max(D);
        
        Max_Look(ii,jj) = Index;
        hd_max_look(ii,jj) = value;
        
    end
end
close(h)

%% Output Saving

if(ch == 1)
    
    f_name_1 = strcat(['Max_Look','.bin']); 
    fileandpath_1=[pathname f_name_1];  
    fid_01 = fopen(fileandpath_1,'wb'); 
    fwrite(fid_01,Max_Look', 'float32');  
    fclose(fid_01);
    
    f_name_2 = strcat(['HD_maxL','.bin']);
    fileandpath_2=[pathname f_name_2];
    fid_02 = fopen(fileandpath_2,'wb');
    fwrite(fid_02,hd_max_look', 'float32');
    fclose(fid_02);
    
end
