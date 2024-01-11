function plotObjLabels5_interactive(N_days,path_main,folder_name,filename,N_measure,idx_Measurement)

    [param,param_plot] = GetIniParam();
    if (param.cont == 1)
        str_cont = 'C';
    else
        str_cont = 'D';
    end

    tau_map_all  = cell(N_days,1);
    chi_map_all  = cell(N_days,1);
    Peak_Img_all = cell(N_days,1);
    Mask         = cell(N_days,1);
    for idx_Day = 1:N_days
        path = strcat(path_main,folder_name(idx_Day,1),'/',filename{idx_Day,N_measure(idx_Measurement)});
        if (path(end) == '/')
            path = path(1:end-1);
        end
        
        filename_saving_temp = split(strip((path)),'/');
        filename_saving = strcat(string(filename_saving_temp(end-1)),string(filename_saving_temp(end)));
        filename_saving_results = strcat('Results/RESULTS_EXP',num2str(param.order),str_cont,'_',filename_saving,'_M',num2str(param.Method),'.mat');
        
        load(strcat('Results/Masks/Mask_',folder_name(idx_Day,1),'_',filename{idx_Day,N_measure(idx_Measurement)},'.mat'),'Label_Mask','N_Obj','Label_Mask_full','Img_Big_Full','-mat');
        load(filename_saving_results,'I_avg','Peak_Img','tau_map','chi_map','HorizontalPixels','VerticalPixels')
        tau_map_all {idx_Day} = tau_map;
        chi_map_all {idx_Day} = chi_map; 
        Peak_Img_all{idx_Day} = Peak_Img;
        Mask{idx_Day} = Label_Mask;
    end
    count = 0;
    hbins_tau = 0.1:0.1:15;
    E = hbins_tau;
    FullMask_temp = double((Label_Mask_full>0));

    meanvalue = zeros(N_days,1);
    idx_mask = 1;
    idx_Label = 1;
    gcf = figure(idx_Measurement);
    gcf.Position = [50 50 1000 800];
    ax1 = axes(gcf,'Position',[0.05 0.32 0.6 0.6]);
    ax2 = copyobj(ax1,gcf);



    imagesc(ax1,Img_Big_Full); %% Plot the image
    clim(ax1,[0.1,0.6]);

    Mask_Ptr =  imagesc(ax2,FullMask_temp, 'AlphaData', FullMask_temp);
    [y_max,x_max] = size(FullMask_temp);
    
    colormap(ax1,'gray')
    Mycolormap  = [[0 0 0];[0.5 0.5 0.5];[1 1 1];[0 1 0]];

    colormap(ax2,Mycolormap)
    clim(ax2,[0,3]);
    set(ax1,'XColor', 'none','YColor','none')
    set(ax2,'XColor', 'none','YColor','none')
    
    ax2.UserData = linkprop([ax1,ax2],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});
    ax2.Visible = 'off';
    
    %%% Zoomed in figure
    ax4 = cell(N_days,1);
    ax5 = cell(N_days,1);
    for j = 1:N_days
        ax4{j} = axes('Position',[0.01+(j-1)*0.3 0 0.28 0.28]);hold on;
        ax5{j} = copyobj(ax4{j},gcf);
        set(ax4{j},'XColor', 'none','YColor','none')
        set(ax5{j},'XColor', 'none','YColor','none')
        ax5{j}.UserData = linkprop([ax4{j},ax5{j}],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});          
        ax5{j}.Visible = 'off';
        set(ax5{j}, 'YDir', 'reverse');
    end
    
    %%% Histogram
    ax3 = axes('Position',[0.72 0.5 0.2 0.28]);hold on;
    hx = xlabel('Lifetime (ns)');
    hy = ylabel('Distribution');

    % Set up the callback function for mouse clicks on the axis
    % set(gcf, 'ButtonDownFcn', @mouseClickCallback);
    set(gcf, 'WindowButtonDownFcn', @mouseClickCallback);
    set(gcf, 'Interruptible', 'on', 'BusyAction', 'cancel');
    % Set up the callback function for mouse movement in the figure
    set(gcf, 'WindowButtonMotionFcn', @mouseMoveCallback);
    plotZoomed(idx_Label);
    PlotHisto(idx_Label);
    
    function mouseClickCallback(~,~)
        % Callback function to track mouse clicks
        coordinates = get(ax2, 'CurrentPoint');  % Get current point in the axis
        x_coord = coordinates(1, 1);
        y_coord = coordinates(1, 2);
        % Perform actions based on the clicked coordinates
        if (x_coord > 0 && x_coord < x_max && y_coord > 0 && y_coord < y_max)
            if (Label_Mask_full(round(y_coord),round(x_coord)) > 0)
                idx_Label = Label_Mask_full(round(y_coord),round(x_coord));
                FullMask_temp (Label_Mask_full ~= idx_Label & Label_Mask_full >  0) = 1;
                FullMask_temp (Label_Mask_full == idx_Label) = 3;
                set(Mask_Ptr ,'CData', FullMask_temp );
                plotZoomed(idx_Label);
                PlotHisto(idx_Label);
            end
        end
    end
    
    function mouseMoveCallback(~, ~)
        % Callback function to track mouse movement within the axis
        coordinates = get(ax2, 'CurrentPoint');  % Get current point in the axis
        x_coord = coordinates(1, 1);
        y_coord = coordinates(1, 2);
        if (x_coord > 0 && x_coord < x_max && y_coord > 0 && y_coord < y_max)
            if (Label_Mask_full(round(y_coord),round(x_coord)) > 0)
                idx_temp = Label_Mask_full(round(y_coord),round(x_coord));
                FullMask_temp ((Label_Mask_full ~= idx_temp) & (Label_Mask_full ~= 0) & FullMask_temp ~= 3) = 1;
                FullMask_temp (Label_Mask_full == idx_temp & FullMask_temp ~= 3) = 2;
                set(Mask_Ptr ,'CData', FullMask_temp );
            end
        end
    end

function plotZoomed(idx_Label)
        
    for idx_mask = 1:N_days
        Mask_Temp = Mask{idx_mask};
        Mask_Temp(Mask_Temp ~= idx_Label) = 0;
        Mask_Temp(Mask_Temp == idx_Label) = 1;
        Points_mask = sum(Mask_Temp(:));
        if (~Points_mask)

            cla(ax4{idx_mask});
            cla(ax5{idx_mask});
            continue;
        end
        stats = regionprops(Mask_Temp, 'Centroid');
        center = stats.Centroid;        %%% Center of the Mask for the current label

        tau_Obj = tau_map_all{idx_mask};
        tau_Obj(~Mask_Temp) = NaN;
        [V,E] = histcounts(tau_Obj,param_plot.hbins_tau);
        W = V/sum(V);
        dE = diff(E(1:2));
        tau_mean = sum((E(1:end-1)+dE/2).*W);
        tau_std = sqrt(sum(W.*(E(1:end-1) + dE/2 - tau_mean).^2)); 
        imagesc(ax4{idx_mask},Peak_Img_all{idx_mask}); %% Plot the image
        imagesc(ax5{idx_mask},tau_Obj, 'AlphaData', Mask_Temp);

        xmin = find(sum(Mask_Temp,1), 1, 'first');
        ymin = find(sum(Mask_Temp,2), 1, 'first');
        xmax = find(sum(Mask_Temp,1), 1, 'last');
        ymax = find(sum(Mask_Temp,2), 1, 'last');
        dx = xmax - xmin;
        dy = ymax - ymin;
        pbaspect(ax4{idx_mask},[xmax-xmin ymax-ymin 1])
        pbaspect(ax5{idx_mask},[xmax-xmin ymax-ymin 1]) 
        xlim(ax4{idx_mask},[max(center(1)-dx,0),min(center(1)+dx,size(Peak_Img_all{idx_mask},2))])
        ylim(ax4{idx_mask},[max(center(2)-dy,0),min(center(2)+dy,size(Peak_Img_all{idx_mask},1))])
        colormap(ax4{idx_mask},'gray')
        colormap(ax5{idx_mask},jet(30))
        clim(ax5{idx_mask},[1,3]);    
        title(strcat('Day ',num2str(idx_mask)))
        if (idx_mask ==  N_days)
            c = colorbar(ax5{idx_mask});
            c.Label.String = 'Lifetime (ns)';
            c.Label.FontSize = 14;    
        end
    end
end


function PlotHisto(idx_Label)
    leg_str  = cell(N_days,1);
    cla(ax3);
    for idx_day = 1:N_days
        Mask_Temp = Mask{idx_day};
        Mask_Temp(Mask_Temp ~= idx_Label) = 0;
        Mask_Temp(Mask_Temp == idx_Label) = 1;

        tau_Obj = tau_map_all{idx_day};
        tau_Obj(~Mask_Temp) = NaN;

        [V,E] = histcounts(tau_Obj,hbins_tau);
        W = V/sum(V);
        meanvalue(idx_day) = sum(E(1:end-1).*W);
        h = plot(ax3,E(1:end-1),W,'LineWidth',2);
        leg_str{idx_day} = strcat('Day ',num2str(idx_day));

        if (idx_day == 1)
            title_str = strcat("Mean: ",num2str(meanvalue(idx_day)));
        else
            title_str = strcat(title_str,"->",num2str(meanvalue(idx_day)));
        end

    end
    leg = legend (ax3,string(leg_str),'location','northwest');
    title(ax3,title_str)
    set(ax3,'XLim',[0,5])

    end
    
end