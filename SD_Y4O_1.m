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

clear all;
clc;
%% Parameter Input

rows=input('Number of Rows: ');
cols=input('Number of Columns: ');

look=input('Number of Looks: ');
wsi=input('Enter Window Size: ');

ang_step=input('Angle Increment Step (1/0.5/0.25): ');

ep = 0.00001;

%% File Input Loop

Final=zeros(cols,rows,9);

for aa=1:1:9
    
    [filename, pathname, filterindex] = uigetfile('*.*', 'Image selection');
    fileandpath=[pathname filename];
    
    disp('   ');
    disp(['File:', fileandpath]);
    
    fid=fopen(fileandpath,'rb');
    
    Final(:,:,aa)=fread(fid,[cols rows], 'float32');
    
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

C11=Final(:,:,1)' + ep;
C22=Final(:,:,2)' + ep;
C33=Final(:,:,3)' + ep;

ReC12=Final(:,:,4)';
ImC12=Final(:,:,5)';

C12=complex(ReC12,ImC12) + ep;
C21=conj(C12);

ReC13=Final(:,:,6)';
ImC13=Final(:,:,7)';

C13=complex(ReC13,ImC13) + ep;
C31=conj(C13);

ReC23=Final(:,:,8)';
ImC23=Final(:,:,9)';

C23=complex(ReC23,ImC23) + ep;
C32=conj(C23);

%%

Real_Rotation=zeros(rows,cols);

Final_Peak=zeros(rows,cols);

Final_Offset=zeros(rows,cols);

%% for window processing

wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopi= rows-inci; % Stop row for window processing
stopj= cols-incj; % Stop column for window processing

%%
t_start = tic;
start_time = datestr(now);
h = waitbar(0,'1','Name','Progress Monitor...');

for ii=starti:stopi
    
    perc = (ii/(stopi))*100;
    waitbar(ii/(stopi),h, sprintf('Start Time:  %s \n%6.2f%% of Process Completed... ', start_time, perc));
    
    for jj=startj:stopj
        
        %%
        Bin_C11=C11(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_C11=mean(Bin_C11(:));
        
        Bin_C22=C22(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_C22=mean(Bin_C22(:));
        
        Bin_C33=C33(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_C33=mean(Bin_C33(:));
        
        %%
        Bin_C12=C12(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_C12=mean(Bin_C12(:));
        
        Bin_C21=C21(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_C21=mean(Bin_C21(:));
        
        %%
        Bin_C13=C13(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_C13=mean(Bin_C13(:));
        
        Bin_C31=C31(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_C31=mean(Bin_C31(:));
        
        %%
        Bin_C23=C23(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_C23=mean(Bin_C23(:));
        
        Bin_C32=C32(ii-inci:ii+inci,jj-incj:jj+incj);
        Bin_C32=mean(Bin_C32(:));
        
        Real_Rotation(ii,jj) = (1./4).*atand((-2.*real(Bin_C23))./(Bin_C33 - Bin_C22));
        
        Hellinger_Dist_22 = zeros(1,max(size(-45:ang_step:45)));
        Hellinger_Dist_33 = zeros(1,max(size(-45:ang_step:45)));
        
        temp = 0;
        
        Hellinger_Dist_T = zeros(1,max(size(-45:ang_step:45)));
        
        
        for Rotation=-45:ang_step:45
            
            C_Unrot= [abs(Bin_C11) Bin_C12 Bin_C13;Bin_C21 abs(Bin_C22) Bin_C23;Bin_C31 Bin_C32 abs(Bin_C33)];
            
            Rot=[1 0 0;0 cosd(2*Rotation) sind(2*Rotation);0 -sind(2*Rotation) cosd(2*Rotation)];
            
            temp=temp+1;
            
            inv_Rot = Rot\eye(size(Rot));
            C_Rot=Rot*C_Unrot*inv_Rot;
            
            C_Rot= [abs(C_Rot(1,1)) C_Rot(1,2) C_Rot(1,3);C_Rot(2,1) abs(C_Rot(2,2)) C_Rot(2,3);C_Rot(3,1) C_Rot(3,2) abs(C_Rot(3,3))];
            
            xx = sqrt((abs(C_Unrot(3,3))*abs(C_Rot(3,3))))/(0.5*abs(C_Unrot(3,3)) + 0.5*abs(C_Rot(3,3)));
            
            yy = sqrt((abs(C_Unrot(2,2))*abs(C_Rot(2,2))))/(0.5*abs(C_Unrot(2,2)) + 0.5*abs(C_Rot(2,2)));
            
            Hellinger_Dist_33(1,temp) = abs(1 - (xx).^(look));
            
            Hellinger_Dist_22(1,temp) = abs(1 - (yy).^(look));
            
            
        end
        
        %% Peak Search
        
        Regional_Peaks=imregionalmax(Hellinger_Dist_33);
        r=length(Hellinger_Dist_33);
        
        %%
        if(Regional_Peaks(1,1)==1) % end peaks have not been removed
            Regional_Peaks(1,1)=0; % consider revising if you face any problem
        end
        
        if(Regional_Peaks(1,r)==1)
            Regional_Peaks(1,r)=0;
        end
        %%
        
        temp=zeros(1,r);
        for ss=1:1:r
            
            if(Regional_Peaks(1,ss)==1)
                temp(1,ss)=ss;
            end
        end
        
        temp=temp(temp~=0);
        
        %% Critical step
        
        Len_vec = length(-45:ang_step:45);
        Inv_step = 1/ang_step;
        Add_step = fix(Len_vec/2)+1;
        
        %Fin_ang = (temp_ang - Add_step)./Inv_step; % logic to retrive the
        %angle
        
        %%
        Peak_Angles = (temp - Add_step)./Inv_step;
        P_A_Length=length(Peak_Angles);
        
        if(P_A_Length == 2)
            
            if(Hellinger_Dist_33(1,Peak_Angles(1,1)*Inv_step + Add_step) > Hellinger_Dist_33(1,Peak_Angles(1,2)*Inv_step + Add_step))
                
                Large_Peak = Peak_Angles(1,1);
                Small_Peak = Peak_Angles(1,2);
                
            else
                
                Large_Peak = Peak_Angles(1,2);
                Small_Peak = Peak_Angles(1,1);
                
            end
        end
        
        if (P_A_Length == 1)
            
            xx=Peak_Angles(1,1);
            Large_Peak=xx;
            Small_Peak=xx;
            
        end
        
        if((Hellinger_Dist_33(1,Large_Peak*Inv_step + Add_step) - Hellinger_Dist_22(1,Large_Peak*Inv_step + Add_step)) > 0)
            
            Final_Peak(ii,jj)=Large_Peak;
            Final_Offset(ii,jj) = (Hellinger_Dist_33(1,Large_Peak*Inv_step + Add_step) - Hellinger_Dist_22(1,Large_Peak*Inv_step + Add_step));
            
        else
            
            Final_Peak(ii,jj)=Small_Peak;
            Final_Offset(ii,jj) = (Hellinger_Dist_33(1,Small_Peak*Inv_step + Add_step) - Hellinger_Dist_22(1,Small_Peak*Inv_step + Add_step));
            
        end
        
        %         if (Final_Peak(ii,jj) > 22.5)
        %             ftemp1 = Final_Peak(ii,jj);
        %             Final_Peak(ii,jj) = ftemp1 - 45;
        %
        %         elseif (Final_Peak(ii,jj) < -22.5)
        %             ftemp2 = Final_Peak(ii,jj);
        %             Final_Peak(ii,jj) = ftemp2 + 45;
        %         end
        
        
    end
    
end
close(h);

% Saving Parameters

ch = input('Do you want to save the parameters: ');

if(ch == 1)
    
    [filename, pathname, filterindex] = uigetfile('*.*', 'Folder to Save');
    
    f_name_1 = strcat(['THEETA_Real_Rot','.bin']);
    f_name_2 = strcat(['THEETA_HD','.bin']);
    f_name_3 = strcat(['HD_T33_T22','.bin']);
    
    % path
    
    fileandpath_1=[pathname f_name_1];
    fileandpath_2=[pathname f_name_2];
    fileandpath_3=[pathname f_name_3];
    
    % open an image to write the raw SAR data
    
    fid_01 = fopen(fileandpath_1,'wb');
    fid_02 = fopen(fileandpath_2,'wb');
    fid_03 = fopen(fileandpath_3,'wb');
    
    % write the raw SAR data into the opened file
    
    fwrite(fid_01,Real_Rotation', 'float32');
    fwrite(fid_02,Final_Peak', 'float32');
    fwrite(fid_03,Final_Offset', 'float32');
    
    % close the file
    
    fclose(fid_01);
    fclose(fid_02);
    fclose(fid_03);
    
end


