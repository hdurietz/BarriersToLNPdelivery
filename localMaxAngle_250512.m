% Input data is contained in separate 130x130 pixel tif images channels
% designated C1-C3 present in folder structure
% 'Users/JaneDoe/MATLAB/DATA/exps'


% Get a list of all endosome folders
files = dir(['/Users/JaneDoe/MATLAB/DATA/exps'])
files(1:2) = [] %remove the first two rows
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir]
% Extract only those that are directories.
subFolders = files(dirFlags)
% struct to table
subFolders = struct2table(subFolders)
endosome = subFolders.name
endosome_nr = size(endosome,1) %endosome_nr is the number of rows (folders)

%create one excel file for each endosome that is saved in each endosome
%folder
result_List_C1 = nan(endosome_nr,2)
result_List_C2 = nan(endosome_nr,2)


for g = 1:endosome_nr
    
    
tempEndosome = endosome{3,1}

cd(['Users/JaneDoe/MATLAB/DATA/exps' tempEndosome ]) 


%read image files
fileList_C1 = dir('C1-*.tif')
name_C1 = fileList_C1.name
C1 = imread(name_C1);

fileList_C2 = dir('C2-*.tif')
name_C2 = fileList_C2.name
C2 = imread(name_C2);

fileList_C3 = dir('C3-*.tif')
name_C3 = fileList_C3.name
C3 = imread(name_C3);



    close('all')
    figure;
    imshow(C1, [30, 400])
    figure;
    imshow(C2, [30, 400])
    figure;
    imshow(C3, [30, 400])

    %create binary image segmentation

    bg = prctile(C3, 5, 'all');
    pc95signal = prctile(C3, 95, 'all');
    threshold = (pc95signal-bg)/3 + bg; %can modulate this one depending on signal 
    mask = C3 > threshold;
    figure
    imshow(mask)

    %find center of mass
    s = regionprops(mask,'centroid');
    centerMass = s.Centroid;
    hold on
    plot(centerMass(1),centerMass(2),'b*')
    hold off

    %find 1, 2 or 3 local maxima in C2 (LNP) channel
    maxima_C2 = imregionalmax(C2, 8);
    masked_maxima_C2 = maxima_C2 .* mask;
    figure
    imshow(masked_maxima_C2)
    local_maxes_C2 = C2 .* uint16(masked_maxima_C2); %find local maxes
    C2max = [];
    i=0  %create a loop to identify up to 3 local maxima
    while i<3 && max(local_maxes_C2, [], 'all') > bg + 10
        [M,I] = max(local_maxes_C2(:));
        [I_row, I_col] = ind2sub(size(local_maxes_C2),I);
        C2max = [C2max ;[I_row, I_col]];
        local_maxes_C2(I_row-5:I_row+5, I_col-5:I_col+5) = 0; % remove surrounding local maxima to max #1
        i = i+1;
    end

    

    %find 1,2 or 3 local maxima in C1 (CHMP_galectin) channel
    maxima_C1 = imregionalmax(C1, 8);
    masked_maxima_C1 = maxima_C1 .* mask;
    figure
    imshow(masked_maxima_C1)
    local_maxes_C1 = C1 .* uint16(masked_maxima_C1);
    n_spots = size(C2max, 1);
    
    C1max = [];
    i=0  %create a loop to identify up to 3 local maxima
    while i<3 && max(local_maxes_C1, [], 'all') > bg + 5
        [M,I] = max(local_maxes_C1(:));
        [I_row, I_col] = ind2sub(size(local_maxes_C1),I);
        C1max = [C1max ;[I_row, I_col]];
        local_maxes_C1(I_row-5:I_row+5, I_col-5:I_col+5) = 0; % remove surrounding local maxima to max #1
        i = i+1;
    end
    

    %create vectors
    vector_C2 = C2max - [centerMass(2), centerMass(1)];
    vector_C2 = [vector_C2(:,2), -vector_C2(:,1)]; %convert matrix indices to cartesian x, y coordinates
    vector_C1 = C1max - [centerMass(2), centerMass(1)];
    vector_C1 = [vector_C1(:,2), -vector_C1(:,1)]; %convert matrix indices to cartesian x, y coordinates

    %calculate the angle of the vector
    polar_C2 = [9*sqrt(vector_C2(:,1).^2 + vector_C2(:,2).^2), rad2deg(atan2(vector_C2(:,2),vector_C2(:,1)))];
    polar_C1 = [9*sqrt(vector_C1(:,1).^2 + vector_C1(:,2).^2), rad2deg(atan2(vector_C1(:,2),vector_C1(:,1)))];
    
    %calculate number of local maxima for exporting list of data
    polar_C2_size = size(polar_C2,1)
    polar_C1_size = size(polar_C1,1)
    
    %create 3x4 matrix and fill it with NaN
    polar_C1_C2 = NaN(3,4);
    
    %export data
    polar_C1_C2(1:polar_C2_size,1:2) = polar_C2; %LNP in column 1 and 2
    polar_C1_C2(1:polar_C1_size,3:4) = polar_C1; %chmp in column 3 and 4
     
    xlswrite('result_List',polar_C1_C2)
    imwrite(mask,'endosomeMask.png')
    imwrite(masked_maxima_C2,'masked_maxima_LNP.png')
    imwrite(masked_maxima_C1,'masked_maxima_chmp.png')
    
end
