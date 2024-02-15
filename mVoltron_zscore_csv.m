%A script to auto identify ROIs based on brigthness and plot the resulting activity

% Loading .ND2 files to MATLAB
mex -setup
loadlibrary('Nd2ReadSdk','Nd2ReadSdk.h')
% interval to be plotted
interval = 5;

folder = uigetdir;

files = dir(fullfile(folder,'*.nd2'));


for m = 1:numel(files)

    filename = files(m).name; %enter whatever the file name is

    filepath = fullfile(folder,filename);
    %get the image info
    [ImageInfo] = ND2Info(filepath);
    %Read the ND2
    [Image] = ND2ReadSingle(filepath);

    %motion correction
    output_filename =  [filename '_reg.mat'];
    options = NoRMCorreSetParms('d1',size(Image,1),'d2',size(Image,2),'init_batch',50,'make_avi',1);
    [M_final,shifts,template,options,col_shift] = normcorre(Image,options);

    %save the motion corrected file
    save(fullfile(folder,output_filename),'M_final');

    %run some simple k-means to find active regions
    L = imsegkmeans((mean(M_final,3)),2); %divide it into 2 clust
    [val, idx] = max(mean(M_final,3),[],'all','linear'); %pick the cluster with the brightest pixel in it

    %extract clusters
    t = arrayfun(@(x) squeeze(mean(bsxfun(@times,M_final,(L==x)),[1 2])),1:2,'UniformOutput',false);

    roi_ind = L(idx);
    %create a mask based on picked ROI
    mask = ones(size(L));
    mask(L==roi_ind) = 2;

    %grab the time trace
    neuron = t{roi_ind};

    clf
    mean_im = mean(M_final,3);
    fin_im = mat2gray(mean_im,[0 10000]);
    rgbImage = cat(3, fin_im, fin_im, fin_im);
    rgbImage(:,:,1) = rgbImage(:,:,1).*mask;
    imshow(rgbImage);
    exportgraphics(gcf,fullfile(folder,strcat(filename,'_im.png')),'Resolution',300);
    save(fullfile(folder,strcat(filename,'.mat')),'neuron','mask')


    nneuron = neuron*-1;
    description = strsplit(ImageInfo.description,' ');
    exposure_time_idx = find(cellfun(@(x) strcmp('Exposure:',x),description))+1;
    exposure_in_ms = 1.17755;
    if exposure_in_ms >200
        exposure_in_ms = exposure_in_ms/1000;
    end
    length_in_seconds = ImageInfo.Attributes.sequenceCount*exposure_in_ms/1000;

    % Photobleaching correction
    time = (0:length(nneuron)-1)' * exposure_in_ms / 1000;

    % Check for NaN or Inf values and replace them with the median value of the neuron signal
    nneuron_no_nan_inf = nneuron;
    nneuron_no_nan_inf(isnan(nneuron) | isinf(nneuron)) = median(nneuron, 'omitnan');

    % Provide initial starting points for the exponential decay model
    % startPoints = [max(nneuron_no_nan_inf), -0.1];

    % Fit the model
    [fitresult, ~] = fit(time, nneuron_no_nan_inf, 'poly9');

    % Subtract the decay curve from the original signal
    photobleaching_corrected_neuron = nneuron - fitresult(time);

    % Make csv file name
    csv_filename = [filename '_csv.csv'];

    % Export photobleaching corrected trace as .csv file
    writematrix(photobleaching_corrected_neuron, csv_filename);

    % turn fluorescent signal into z score
    z_corr_neur = zscore(photobleaching_corrected_neuron);

    % Further filter the signal
    [b, a] = butter(2, 0.5);
    nfneuron = filtfilt(b, a, double(z_corr_neur));

    % Apply a bandpass filter to isolate 2-10Hz oscillations
    [b_bp, a_bp] = butter(2, [1, 8] / (0.5 * (1 / (exposure_in_ms / 1000))));
    oscillations = filtfilt(b_bp, a_bp, double(nfneuron));

    % Calculate the number of rows needed to cover the entire recording in 5s increments
    rows = ceil(length_in_seconds / 5);

% Plotting with spike indicators (raster marks) and 2-6Hz oscillations
clf
for i = 0:(rows - 1)
    ax = subplot(rows, 1, i + 1);
    
    % Define the start and end points for each 5s increment
    startPoint = round(i * 5 * 1000 / exposure_in_ms) + 1;
    endPoint = round((i + 1) * 5 * 1000 / exposure_in_ms);
    if endPoint > length(nfneuron)
        endPoint = length(nfneuron);
    end
    
    % Plot the original signal
    plot(startPoint:endPoint, nfneuron(startPoint:endPoint))
    
    xticks(startPoint:1000 / exposure_in_ms * 5:endPoint)
    xticklabels(0:5:(endPoint - startPoint) * exposure_in_ms / 1000)
    xlim([startPoint, endPoint])
    set(gca, 'TickDir', 'out')
    ylabel('Fluorescence a.u.')
    
    if i == (rows - 1)
        xlabel(['Time (s)'])
        ax.Units = 'normalized';
        ax.Position = ax.Position + [0, -0.05, 0, 0.05]; % Adjust the position of the bottom subplot
    end
    box off
end

title(strcat('5s increments ', filename), 'Interpreter', 'none')
export_name = strcat(filename, '_5s_increments.pdf');
exportgraphics(gcf, fullfile(folder, export_name), 'Resolution', 300, 'ContentType', 'vector');



    clf
    plot(nfneuron)
    pbaspect([8,2,1])
    title(strcat('Single_plot',filename),'Interpreter','none')
    xticks(0:1000/exposure_in_ms*5:length(nfneuron))
    xticklabels(0:5:length(nfneuron)*exposure_in_ms/1000)
    set(gca,'TickDir','out')
    ylabel('Fluorescence a.u.')
    xlabel('Time (s)')
    export_name = strcat(filename,'_single_plot.pdf');
    box off
    exportgraphics(gcf,fullfile(folder,export_name),'Resolution',300,'ContentType','vector');
end
%end

