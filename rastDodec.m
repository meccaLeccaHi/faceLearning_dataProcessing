function rastDodec(paths,monkey,suffix)
%   loop through relevant files in 'paths.mas'
%
% Usage: rastDodec(paths,monkey,suffix)
% MONKEY is the monkey's name
% PATHS project path structure containing at minimum paths.rare field
% SUFFIX is a string appended to the end of each filename
%
% last modified 6-17-16
% apj

% % housekeeping
warning('off','MATLAB:LargeImage')
% warning('off','MATLAB:opengl:switchToSoftwareUnixNotSupported')
% opengl software

% setup
visTog                                          = 'off';
% lw              = 1;          % LineWidth
% fsz             = 8;          % axes fontsize
% fileFsz         = 13;         % File label fontsize
% dateFsz         = 12;         % Date label fontsize
% titleFsz        = 8;          % Raster title fontsize
% % msz            = 6;          % MarkerSize
% plLw            = 2;          % LineWidth for sdf's
Plotwin                                         = [-300 600]; % plot display size (msec)
Rspwin                                          = [100 300]; % time window for evoked rate
prs                                             = 150;       % print resolution
fullscr                                         = get(0,'ScreenSize'); % get screen size

[x,y,z]                                         = sphere;

dodecStr                                        = fullfile(paths.results,'dodecs','quadPlots');
% faceStim = 8:8:192;
% normStim = [193 194];

% get file list
exp_code                                        = [monkey(1) paths.code];
dirList                                         = dir(fullfile(paths.mas,[exp_code '*' suffix.load '.mat']));


% load all image files for disp (stimscreen 1 & 2)
load(fullfile(paths.stim,'imgFaceLearningSet1OnWhite.mat'));
stim_images                                     = img;
% stim_images = load([paths.stim 'imgFaceLearningSet1.mat']);

species                                         = {'human' 'monkey'};

% loop through file list
for i = 1:length(dirList);
    %     i = 4
    
    % load data file
    filename                                    = fullfile(paths.mas,dirList(i).name);
    load(filename)
    disp(['Loaded: ' dirList(i).name])
    unpack
    
    dates                                       = datenum(dat.h.date,'dd-mmm-yy');
    morphs                                      = dat.h.stim.morphSteps;
    morphList                                   = [0 morphs(5:end)];
    
    % get background rate for each neuron
    trialBackGround                             = [];
    for sc = 1:length(dat.s)
        backGroundSpikeCount                    = length(find(dat.s(sc,:)<-100));
        trialBackGround                         = [trialBackGround; backGroundSpikeCount];
    end
    
    ave_bResp                                   = mean(trialBackGround)*(1000/diff(Rspwin));
    
    % loop through face types (human v. monkey)
    for typeNum = 1:2
        %         typeNum = 2
        
        % loop through groups of 4 faces
        for quadGroup = 1:3
            
            %% get response range
            
            % set face numbers
            tempFaceList                        = 100*(typeNum-1)+(1:4)+4*(quadGroup-1);
            %             tempFaceList = [104 102 106 109]
            
            % compare responses to all faces of each species
            maxResp                             = -99;
            minResp                             = 99;
            for ii = 1:length(tempFaceList)
                fNum                            = [100*(typeNum-1) tempFaceList(ii).*ones(1,length(morphs(5:end)))];
                for iii = 1:length(fNum)
                    stimDatStep = select_trials(dat,TYPE,typeNum,FACE,fNum(iii),...
                        STEP,morphList(iii));
                    trialSpikes = [];
                    for sc = 1:length(stimDatStep.s(:,1))
                        trialSpikeCount         = length(find(stimDatStep.s(sc,:)>...
                            Rspwin(1)&stimDatStep.s(sc,:)<Rspwin(2)));
                        trialSpikes             = [trialSpikes; trialSpikeCount];
                    end
                    eResp                       = mean(trialSpikes)*(1000/diff(Rspwin));
                    tResp                       = eResp-ave_bResp;
                    if tResp>maxResp
                        maxResp                 = tResp;
                    end
                    if tResp<minResp
                        minResp                 = tResp;
                    end
                end
            end
            
            %% set color axis according to response range
            ncolors                             = 32;
            cc                                  = flipud(hot(ncolors));
            cutoffs                             = nan(ncolors,2);
            level                               = (maxResp-minResp)/ncolors;
            % loop through each color in range
            for u = 1:ncolors
                cutoffs(u,1)                    = str2double(sprintf('%4.2f',maxResp-((maxResp-minResp)*u/ncolors)));
                cutoffs(u,2)                    = str2double(sprintf('%4.2f',maxResp-((maxResp-minResp)*u/ncolors)+level));
            end
            cutoffs(1,2)                        = str2double(sprintf('%4.2f',maxResp+abs(.1*maxResp)));
            cutoffs(end,1)                      = str2double(sprintf('%4.2f',minResp-abs(.1*maxResp)));
            
            if maxResp>=1
                
                %% plot dodec
                figure(27);clf
                figPos                          = fullscr.*[1 1 .9 .9];
                if regexp(visTog,'off')
                    figPos(1)                   = figPos(1)-1000;
                end
                set(27,'Position',figPos,'Visible',visTog);
                
                caxis([minResp maxResp]);
                axis([-1 1 -1 1]); hold on
                aH                              = gca;
                axis off
                mainAxPos                       = get(aH,'Position');
                
                % get max/min values (for this group of faces) to normalize sdf's
                tempMaxVec                      = [];
                tempMinVec                      = [];
                for f = 1:length(tempFaceList)
                    fNum                        = [100*(typeNum-1) tempFaceList(f).*ones(1,length(morphs(5:end)))];
                    for g = 1:length(fNum)
                        stepDatRast             = select_trials(dat,TYPE,typeNum,FACE,fNum(g),STEP,morphList(g));
                        tempSdf                 = get_sdf(stepDatRast,1,Plotwin);
                        tempMaxVec              = [tempMaxVec; max(tempSdf)];
                        tempMinVec              = [tempMinVec; min(tempSdf)];
                    end
                end
                sdfMax                          = max(tempMaxVec);
                sdfMin                          = min(tempMinVec);
                
                for f = 1:length(tempFaceList)
                    
                    fNum                        = [100*(typeNum-1) tempFaceList(f).*ones(1,length(morphs(5:end)))];
                    
                    axWdth                      = .05;
                    pH2                         = axWdth/2;
                    switch f
                        case 1
                            flipY                   = 1;
                            placeHolder             = 0;
                            skootch                 = .015;
                            xFlip                   = 1;
                        case 2
                            flipY                   = -1;
                            placeHolder             = .005;
                            skootch                 = 0;
                            xFlip                   = 1;
                        case 3
                            flipY                   = -1;
                            placeHolder             = .005;
                            skootch                 = 0;
                            xFlip                   = -1;
                        case 4
                            flipY                   = 1;
                            placeHolder             = 0;
                            skootch                 = .015;
                            xFlip                   = -1;
                    end
                    
                    for g = 1:length(fNum)
                        
                        set(27,'CurrentAxes',aH)
                        
                        stepDat                 = select_trials(dat,TYPE,typeNum,FACE,fNum(g),STEP,morphList(g));
                        imNum                   = unique(stepDat.c(:,STIM)); % stimulus image index
                        stepSize                = .3/length(fNum);
                        
                        trialSpikes = [];
                        for sc = 1:length(stepDat.s(:,1))
                            trialSpikeCount = length(find(stepDat.s(sc,:)>Rspwin(1)&stepDat.s(sc,:)<Rspwin(2)));
                            trialSpikes         = [trialSpikes; trialSpikeCount];
                        end
                        eResp                   = mean(trialSpikes)*(1000/diff(Rspwin));
                        tResp                   = eResp-ave_bResp; % sprintf('%4.2f',eResp-ave_bResp)
                        
                        col                     = find(tResp>=cutoffs(:,1)&tResp<=cutoffs(:,2),1,'first');
                        
                        tempIm                  = stim_images{imNum(end)};
                        
                        %                         [imHt,imWidth,~]    = size(tempIm);
                        %                         axHt                = (imHt/imWidth)*axWdth;
                        
                        %         [az,el] = view
                        %         view(az,el);
                        
                        % plot sphere
                        tempPos                 = [mainAxPos(1)+mainAxPos(3)/2+(g-1)*stepSize*xFlip...
                            mainAxPos(2)+mainAxPos(4)/2.5+(g-1)*stepSize/2*flipY...
                            mainAxPos(3)/10 mainAxPos(4)/10];
                        axes('Position',tempPos);
                        axis([0 1 0 1 0 1])
                        b                       = surf(x,y,z,'Facecolor',cc(col,:),'LineStyle','none');
                        light % lightangle(-45,30);
                        set(b,'DiffuseStrength',.9,'AmbientStrength',.6 ,...
                            'EdgeLighting','gouraud','FaceLighting','gouraud');
                        
                        %             r = .5;x = 0.5;y = 0.5;
                        %             circles(x,y,r,'EdgeColor',cc(col,:),'Color',cc(col,:))
                        axis square off
                        
                        if g>=2||f==1
                            
                            % plot face image
                            imShrink            = (typeNum-1)*.02;
                            tempPos             = [mainAxPos(1)+mainAxPos(3)/2+(g-1)*stepSize*xFlip...
                                (4*placeHolder*(typeNum-1))+mainAxPos(2)+mainAxPos(4)/2.5+mainAxPos(4)/10*flipY+(g-1)*stepSize/2*flipY ...
                                mainAxPos(3)/10 mainAxPos(4)/10-imShrink];
                            axes('Position',tempPos);
                            tempRS              = imresize(tempIm,.5);
                            image(tempRS)
                            axis image off
                            
                            % plot raster
                            if length(dates)~=1
                                dates           = dates(1:3);
                            end
                            stepDatRast         = select_trials(dat,TYPE,typeNum,FACE,fNum(g),STEP,morphList(g),DATE,dates);
                            tempPos             = [pH2/2+mainAxPos(1)+mainAxPos(3)/2+(g-1)*stepSize*xFlip...
                                placeHolder-7*placeHolder+pH2/2+mainAxPos(2)+mainAxPos(4)/2.5+mainAxPos(4)/5*flipY+(g-1)*stepSize/2*flipY-imShrink*flipY ...
                                axWdth axWdth];
                            tempPos(2)          = tempPos(2)-skootch;
                            axes('Position',tempPos);
                            plot_rasters_dark(stepDatRast.s,Plotwin,1.5);
                            %                             if f==1|f==4
                            %                                 tempYLim = ylim;
                            %                                 ylim([min(tempYLim) 30])
                            %                                 if g == length(fNum)
                            %                                     line([0 300],[10 10],'LineWidth',10,'Color','k')
                            %                                 end
                            %                             else
                            tempYLim            = ylim;
                            ylim([min(tempYLim)-30 0])
                            if g == length(fNum)
                                line([30 270],[min(tempYLim)-5 min(tempYLim)-5],'LineWidth',6,'Color','k')
                                if f==1
                                    text(-35,min(tempYLim)-21,'0','FontWeight','bold','Fontsize',8)
                                    text(170,min(tempYLim)-21,'300','FontWeight','bold','Fontsize',8)
                                    text(75,min(tempYLim)-32,'ms','Fontsize',7.5)
                                    %                                 elseif f==2
                                    %                                     notes_pos   = get(gca,'Position');
                                end
                            end
                            %                             end
                            axis off
                            
                            % plot sdf
                            tempSdf             = get_sdf(stepDat,1,Plotwin);
                            tempPos(2)          = tempPos(2)+.0505;
                            axes('Position',[tempPos(1:3) .0475]);
                            area(1:length(tempSdf),tempSdf,'EdgeColor',[.7 .7 .7],'FaceColor',[.7 .7 .7],...
                                'LineWidth',1); hold on
                            %                             plot(tempSdf,'LineWidth',2)
                            xlim([1 length(tempSdf)])
                            ylim([sdfMin sdfMax])
                            axis off
                        end
                    end
                end
                
                set(27,'CurrentAxes',aH,'Color','w','Visible',visTog)
                xDrop                           = .1;
                yDrop                           = .08;
                line([-.65 .65]+xDrop,[-.33 .295]-yDrop,'Color',[.9 .9 .9],'LineWidth',9,'LineStyle','-')
                line([-.65 .65]+xDrop,[.295 -.33]-yDrop,'Color',[.9 .9 .9],'LineWidth',9,'LineStyle','-')
                
                colPos                          = [.5 .2 .11 .015];
                ha                              = axes('Position',colPos);
                set(ha,'XTick',[],'YTick',[])
                colormap(hot);
                hcb                             = colorbar('Location','North','Position',colPos);
                tempLabel                       = str2num(num2str(round([minResp(1) maxResp].*100)/100,3));
                set(get(hcb,'XLabel'),'String','{\itHz}','FontWeight','bold');
                set(hcb,'XTick',get(hcb,'xlim'),'XTickLabel',tempLabel,...
                    'FontWeight','bold')
                % 'YTick',[];
                
                filename                        = [dirList(i).name(1:end-4)...
                    '_' species{typeNum} 'faceGroup' num2str(quadGroup) suffix.save '.png'];
                savename                        = fullfile(dodecStr,monkey,filename);
                %                 savename = [dodecStr 'quad' dirList(i).name(1:end-4) 'type' num2str(typeNum) 'Parade.png'];
                export_fig(27,savename,['-r' num2str(prs)]); % ,'-nocrop'
                disp(['Saved figure: ' filename])
                clf
            else
                disp(['Skipped: Max firing rate for face group #' num2str(quadGroup)...
                    ' is too low (' num2str(round(maxResp*100)/100) 'Hz)']);
            end
        end
    end
end
end