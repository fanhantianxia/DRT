function wb_plot_networks(EEG_results,outputdir,filename)
% plot results of networks calculated by wb_calculate_EEGnetwork
% Input£º
%    EEG_results:structure EEG_results calculated by wb_calculate_EEGnetwork
%    outputdir: save dir of figures
%    filename: data file name (subject)
% -------------------------------------------------------------------------
% Written by Li Dong (UESTC,Lidong@uestc.edu.cn)
% $ 2018.5.29 £¨NOT completed£©
% -------------------------------------------------------------------------
if ~isempty(EEG_results) && isfield(EEG_results,'type')
    if isequal(EEG_results.type,'network')
        if ~isempty(EEG_results.parameter.chanlocs)
            chanlocs = EEG_results.parameter.chanlocs(EEG_results.parameter.selechanns);
            % Power_mean
            if ~isempty(EEG_results.M_mean)
                try
                    h1 = figure('visible','off');
                    set(h1, 'position', get(0,'ScreenSize'));
                    DIM = size(EEG_results.M_mean);
                    if length(DIM) == 2
                    elseif length(DIM) == 3
                        Nf = DIM(3);
                        for j1 = 1:Nf
                            if any(~isnan(EEG_results.M_mean(:,:,j1)))
                                subplot(ceil(Nf/4),4,j1);
                                topoplot(EEG_results.M_mean(:,:,j1),chanlocs,'shrink','off','maplimits','maxmin');
                                axis equal;
                                colorbar;colormap(jet);
                                title([EEG_results.parameter.bandName{j1},...
                                    ':',num2str(EEG_results.parameter.bandLimit(j1,1)),...
                                    '-',num2str(EEG_results.parameter.bandLimit(j1,2)),'Hz']);
                            end;
                        end;
                    end
                    set(h1,'inverthardcopy','off','PaperPositionMode','auto');
                    print(h1,'-dpng','-r300',fullfile(outputdir,[filename,'-EEG-results-Power_mean.png']));
                catch
                    disp('failed to plot results')
                end;
            end
            
            % Power_relative_mean
            if ~isempty(EEG_results.Power_relative_mean)
                try
                    h2 = figure('visible','off');
                    set(h2, 'position', get(0,'ScreenSize'));
                    Nf = size(EEG_results.Power_relative_mean,2);
                    for j1 = 1:Nf
                        if any(~isnan(EEG_results.Power_relative_mean(:,j1)))
                            subplot(ceil(Nf/4),4,j1);
                            topoplot(EEG_results.Power_relative_mean(:,j1),chanlocs,'shrink','off','maplimits','maxmin');
                            axis equal;
                            colorbar;colormap(jet);
                            title(['Relative ',EEG_results.parameter.bandName{j1},...
                                ':',num2str(EEG_results.parameter.bandLimit(j1,1)),...
                                '-',num2str(EEG_results.parameter.bandLimit(j1,2)),'Hz']);
                        end
                    end;
                    set(h2,'inverthardcopy','off','PaperPositionMode','auto');
                    print(h2,'-dpng','-r300',fullfile(outputdir,[filename,'-EEG-results-Power_relative_mean.png']));
                catch
                    disp('failed to plot results')
                end;
            end;
            
            % mean R1-R6, PAF
            try
                h3 = figure('visible','off');
                set(h3, 'position', get(0,'ScreenSize'));
                
                if any(~isnan(EEG_results.R1_mean))
                    subplot(2,4,1);
                    topoplot(EEG_results.R1_mean,chanlocs,'shrink','off','maplimits','maxmin');
                    axis equal;
                    colorbar;colormap(jet);
                    title('R1 mean');
                end
                
                if any(~isnan(EEG_results.R2_mean))
                    subplot(2,4,2);
                    topoplot(EEG_results.R2_mean,chanlocs,'shrink','off','maplimits','maxmin');
                    axis equal;
                    colorbar;colormap(jet);
                    title('R2 mean');
                end
                
                if any(~isnan(EEG_results.R3_mean))
                    subplot(2,4,3);
                    topoplot(EEG_results.R3_mean,chanlocs,'shrink','off','maplimits','maxmin');
                    axis equal;
                    colorbar;colormap(jet);
                    title('R3 mean');
                end
                
                if any(~isnan(EEG_results.R4_mean))
                    subplot(2,4,4);
                    topoplot(EEG_results.R4_mean,chanlocs,'shrink','off','maplimits','maxmin');
                    axis equal;
                    colorbar;colormap(jet);
                    title('R4 mean');
                end
                
                if any(~isnan(EEG_results.R5_mean))
                    subplot(2,4,5);
                    topoplot(EEG_results.R5_mean,chanlocs,'shrink','off','maplimits','maxmin');
                    axis equal;
                    colorbar;colormap(jet);
                    title('R5 mean');
                end
                
                if any(~isnan(EEG_results.R6_mean))
                    subplot(2,4,6);
                    topoplot(EEG_results.R6_mean,chanlocs,'shrink','off','maplimits','maxmin');
                    axis equal;
                    colorbar;colormap(jet);
                    title('R6 mean');
                end
                
                if any(~isnan(EEG_results.PAF_mean))
                    subplot(2,4,7);
                    topoplot(EEG_results.PAF_mean,chanlocs,'shrink','off','maplimits','maxmin');
                    axis equal;
                    colorbar;colormap(jet);
                    title('PAF mean');
                end
                
                set(h3,'inverthardcopy','off','PaperPositionMode','auto');
                print(h3,'-dpng','-r300',fullfile(outputdir,[filename,'-EEG_results-R1-R6-PAF.png']));
            catch
                disp('failed to plot results')
            end;
        else
            if ~isempty(EEG_results.Power_mean)
                try
                    h1 = figure('visible','off');
                    set(h1, 'position', get(0,'ScreenSize'));
                    Nf = size(EEG_results.Power_mean,2);
                    for j1 = 1:Nf
                        if any(~isnan(EEG_results.Power_mean(:,j1)))
                            subplot(ceil(Nf/4),4,j1);
                            plot(EEG_results.Power_mean(:,j1));
                            xlabel('channels');
                            ylabel('power');
                            title([EEG_results.parameter.bandName{j1},...
                                ':',num2str(EEG_results.parameter.bandLimit(j1,1)),...
                                '-',num2str(EEG_results.parameter.bandLimit(j1,2)),'Hz']);
                        end;
                    end;
                    set(h1,'inverthardcopy','off','PaperPositionMode','auto');
                    print(h1,'-dpng','-r300',fullfile(outputdir,[filename,'-EEG-results-Power_mean.png']));
                catch
                    disp('failed to plot results')
                end;
            end
            
            % Power_relative_mean
            if ~isempty(EEG_results.Power_relative_mean)
                try
                    h2 = figure('visible','off');
                    set(h2, 'position', get(0,'ScreenSize'));
                    Nf = size(EEG_results.Power_relative_mean,2);
                    for j1 = 1:Nf
                        if any(~isnan(EEG_results.Power_relative_mean(:,j1)))
                            subplot(ceil(Nf/4),4,j1);
                            plot(EEG_results.Power_relative_mean(:,j1));
                            xlabel('channels');
                            ylabel('relative power');
                            title(['Relative ',EEG_results.parameter.bandName{j1},...
                                ':',num2str(EEG_results.parameter.bandLimit(j1,1)),...
                                '-',num2str(EEG_results.parameter.bandLimit(j1,2)),'Hz']);
                        end
                    end;
                    set(h2,'inverthardcopy','off','PaperPositionMode','auto');
                    print(h2,'-dpng','-r300',fullfile(outputdir,[filename,'-EEG-results-Power_relative_mean.png']));
                catch
                    disp('failed to plot results')
                end;
            end;
            
            % mean R1-R6, PAF
            try
                h3 = figure('visible','off');
                set(h3, 'position', get(0,'ScreenSize'));
                
                if any(~isnan(EEG_results.R1_mean))
                    subplot(2,4,1);
                    plot(EEG_results.R1_mean);
                    xlabel('channels');
                    ylabel('R1');
                    title('R1 mean');
                end
                
                if any(~isnan(EEG_results.R2_mean))
                    subplot(2,4,2);
                    plot(EEG_results.R2_mean);
                    xlabel('channels');
                    ylabel('R2');
                    title('R2 mean');
                end
                
                if any(~isnan(EEG_results.R3_mean))
                    subplot(2,4,3);
                    plot(EEG_results.R3_mean);
                    xlabel('channels');
                    ylabel('R3');
                    title('R3 mean');
                end
                
                if any(~isnan(EEG_results.R4_mean))
                    subplot(2,4,4);
                    plot(EEG_results.R4_mean);
                    xlabel('channels');
                    ylabel('R4');
                    title('R4 mean');
                end
                
                if any(~isnan(EEG_results.R5_mean))
                    subplot(2,4,5);
                    plot(EEG_results.R5_mean);
                    xlabel('channels');
                    ylabel('R5');
                    title('R5 mean');
                end
                
                if any(~isnan(EEG_results.R6_mean))
                    subplot(2,4,6);
                    plot(EEG_results.R6_mean);
                    xlabel('channels');
                    ylabel('R6');
                    title('R6 mean');
                end
                
                if any(~isnan(EEG_results.PAF_mean))
                    subplot(2,4,7);
                    plot(EEG_results.PAF_mean);
                    xlabel('channels');
                    ylabel('PAF');
                    title('PAF mean');
                end
                
                set(h3,'inverthardcopy','off','PaperPositionMode','auto');
                print(h3,'-dpng','-r300',fullfile(outputdir,[filename,'-EEG_results-R1-R6-PAF.png']));
            catch
                disp('failed to plot results')
            end;
        end
    else
        warning('the type of results is not correct')
    end
else
    disp('EEG_reuslts is empty or the type of results is unknown')
end
