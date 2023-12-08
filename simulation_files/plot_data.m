% tarot test flight data plot

% demux muxed data 
% adjust jitter in time values when close to sample time 
Tt = T;
for i=2:length(T),
    deltat = T(i)- Tt(i-1); 
    if deltat > 0.022 && deltat < 0.027,        % se supero il tempo di camp
        %no missing data BUT sample time a little bit longer -> restore good value
        Tt(i) = Tt(i-1)+0.02;
    elseif deltat < 0.018,
        %no missing data BUT sample time a little bit shorter -> restore good value
        Tt(i) = Tt(i-1)+0.02;
    elseif deltat >= 0.027 
        %missing data -> must find exact value that is multiple of sample time
        count = 0;
        while deltat >= 0.018,
            deltat = deltat-0.02;
            count = count+1;
        end
        Tt(i) = Tt(i-1)+count*0.02;
    else
        %this is a good sample already -> do nothing
    end
end
%debug plot
figure,
plot(diff(T),'bo'), hold on, plot(diff(Tt),'rs')
% demux
T_100 = Tt.*100;
t_100_demux = mod(T_100,4);
figure, plot(t_100_demux,'bo');
%
x_idx  = find(t_100_demux < 1 | t_100_demux >3);
y_idx  = find(t_100_demux > 1 & t_100_demux <3);
%
% now length(x_idx)+length(y_idx) must be = length(T)

% SET FIRMWARE REVISION
%
%1 - firmware version with debug_pos_control only (muxed data in x and y at positions 41:45
%2 - firmware version with debug_pos_control (41:45) debug_h_control and debug_waypoints (35:39)
firmware_version = 2;
%
% and process
if firmware_version == 1,
    %demux debug_pos_control 
    DEBUG_POS_CONTROL_x = DATA(:,41:45);
    DEBUG_POS_CONTROL_x(y_idx,:) = NaN;
    DEBUG_POS_CONTROL_y = DATA(:,41:45);
    DEBUG_POS_CONTROL_y(x_idx,:) = NaN;
    %replace NANs
    for i = 1:length(T),
        if isnan(DEBUG_POS_CONTROL_x(i,1)), if i>1, DEBUG_POS_CONTROL_x(i,:) = DEBUG_POS_CONTROL_x(i-1,:); else DEBUG_POS_CONTROL_x(i,:)=0; end; end;
        if isnan(DEBUG_POS_CONTROL_y(i,1)), if i>1, DEBUG_POS_CONTROL_y(i,:) = DEBUG_POS_CONTROL_y(i-1,:); else DEBUG_POS_CONTROL_y(i,:)=0; end; end;
    end    
    PILOT_COMMANDS = DATA(:,35:38);
    MODE = DATA(:,39);
elseif firmware_version == 2,
    %demux debug_pos_control 
    DEBUG_POS_CONTROL_x = DATA(:,41:45); DEBUG_POS_CONTROL_x(y_idx,:) = NaN;
    DEBUG_POS_CONTROL_y = DATA(:,41:45); DEBUG_POS_CONTROL_y(x_idx,:) = NaN;
    %demux pilot commands, mode and debug_waypoints
    PILOT_COMMANDS = DATA(:,35:38); PILOT_COMMANDS(x_idx,:) = NaN; 
    MODE = DATA(:,39); MODE(x_idx) = NaN;
    DEBUG_WAYP = DATA(:,35:39); DEBUG_WAYP(y_idx,:) = NaN;
    %replace NANs
    for i = 1:length(T),
        if isnan(DEBUG_POS_CONTROL_x(i,1)), if i>1, DEBUG_POS_CONTROL_x(i,:) = DEBUG_POS_CONTROL_x(i-1,:); else DEBUG_POS_CONTROL_x(i,:)=0; end; end;
        if isnan(DEBUG_POS_CONTROL_y(i,1)), if i>1, DEBUG_POS_CONTROL_y(i,:) = DEBUG_POS_CONTROL_y(i-1,:); else DEBUG_POS_CONTROL_y(i,:)=0; end; end;
        if isnan(PILOT_COMMANDS(i,1)), if i>1, PILOT_COMMANDS(i,:) = PILOT_COMMANDS(i-1,:); else PILOT_COMMANDS(i,:)=0; end; end;
        if isnan(MODE(i,1)), if i>1, MODE(i,:) = MODE(i-1,:); else MODE(i,:)=0; end; end;
        if isnan(DEBUG_WAYP(i,1)), if i>1, DEBUG_WAYP(i,:) = DEBUG_WAYP(i-1,:); else DEBUG_WAYP(i,:)=0; end; end;
    end    
end


%% plot time data
%
% here a+we assume to have T and DATA 
% DATA = DATA_rerun(1:length(T),:);
figure, plot(T, DATA(:,2:4).*180/pi), legend('roll','pitch','yaw');

%inertial data
figure, plot(T, DATA(:,5:7).*180/pi), legend('gyro x','gyro y','gyro z');
figure, plot(T, DATA(:,8:10)), legend('acc x','acc y','acc z');
figure, plot(T, DATA(:,11:13)), legend('mag x','mag y','mag z');
figure, plot(T, sqrt(DATA(:,11).^2+DATA(:,12).^2+DATA(:,13).^2)), legend('mag norm');
%for mag calib
if 0,
    figure, plot(DATA(:,11:13));
    seg = 3000:10300;
    clear MNAV;
    MNAV(:,7:9) = double(DATA(seg,11:13));
end

%ultrasound and baro 
figure, plot(T, DATA(:,14)), legend('ultrasound');
figure, plot(T, DATA(:,15)), legend('h baro filtered');

%navigation data (KF)
figure, plot(T, DATA(:,[1 4 7]+15)), legend('x','y','z');
figure, plot(T, DATA(:,[1 4 7]+15+1)), legend('vx','vy','vz');
figure, plot(T, DATA(:,[1 4 7]+15+2)), legend('ax','ay','az');
figure, plot(T, DATA(:,34)), legend('KF flags');
%compute gps delta T
tgps = T(find(DATA(:,34)==111 | DATA(:,34)==101));
figure, plot(diff(tgps));
disp(['Average GPS fix frequency = ', num2str(1/mean(diff(tgps))), ' Hz']);

% DATI ESPRESSI IN NED (è rimossa la gravità)
%lmeasurement data 
figure, plot(T, DATA(:,[1 4 7]+15+9)), legend('x','y','z'); axis([0 max(T) -150 150]), grid on;
figure, plot(T, DATA(:,[1 4 7]+15+1+9)), legend('vx','vy','vz');
figure, plot(T, DATA(:,[1 4 7]+15+2+9)), legend('ax','ay','az');
figure, plot(T, DATA(:,[1]+15+2+9)), legend('ax');
figure, plot(T, DATA(:,[4]+15+2+9)), legend('ay');
figure, plot(T, DATA(:,[7]+15+2+9)), legend('az');

%altitude
figure, plot(T, DATA(:,15),'r'), hold on;
plot(T, -DATA(:,[7]+15),'b'), 
plot(T, -DATA(:,[7]+15+9),'m'), 
if exist('T_'),
    plot(T_, -DATA_(:,15+7),'g');
    legend('h baro filtered','h kf','h baro','h baro filtered RECOMPUTED');
else
    legend('h kf filtered','h kf','h baro');
end
hold off;
grid on;

% kf position filtering
figure, plot(T, DATA(:,[1]+15),'r','linewidth',4), hold on, plot(T, DATA(:,[4]+15),'g','linewidth',4), plot(T, DATA(:,[7]+15),'b','linewidth',4), 
plot(T, DATA(:,[1]+15+9),'y','linewidth',2), plot(T, DATA(:,[4]+15+9),'m','linewidth',2), plot(T, DATA(:,[7]+15+9),'c','linewidth',2), 
axis([0 max(T) -150 150]), grid on, hold off
legend('estim x','estim y','estim z','gps x','gps y', 'gps z');


%% plot map data

%map view 2D (x=E, y=N)
figure, plot(DATA(:,[4]+15),DATA(:,[1]+15),'r'), hold on 
plot(DATA(:,[4]+15+9),DATA(:,[1]+15+9),'b') 
axis([-150 150 -150 150]), grid on, hold off
legend('estim xy','gps xy');

%map view 2D (x=E, y=N) with coloring based on mode
% divide in sections where mode changes
mode_changes_idx = [1; find(abs(diff(MODE))>0.5); length(T)];
colors = {'k','b','r','m','gx-'};% 
figure;
for i = 1:(length(mode_changes_idx)-1),
    idx_seg = (mode_changes_idx(i)+1):(mode_changes_idx(i+1));
    seg_color = DATA(mode_changes_idx(i)+1,39)+1;
    plot(DATA(idx_seg,[4]+15),DATA(idx_seg,[1]+15),colors{seg_color}), hold on 
end;
hold off
axis equal 
grid on



%map view 3D (x=E, y=N, Z= D)
figure, plot3(DATA(:,[4]+15),DATA(:,[1]+15),-DATA(:,[7]+15),'r'), hold on 
plot3(DATA(:,[4]+15+9),DATA(:,[1]+15+9),-DATA(:,[7]+15+9),'b'), hold off
axis([-150 150 -150 150 -100 50]), grid on, hold off
legend('estim xy','gps xy');

%map view 2D (x=E, y=N) with velocity
figure, plot(DATA(:,[4]+15),DATA(:,[1]+15),'r'), hold on 
plot(DATA(:,[4]+15+9),DATA(:,[1]+15+9),'b') 
%plot vel 
dsf = 50; %down sample factor
quiver(DATA(1:dsf:end,[4]+15),DATA(1:dsf:end,[1]+15),...
    DATA(1:dsf:end,[4]+15+1),DATA(1:dsf:end,[1]+15+1),'g');
quiver(DATA(1:dsf:end,[4]+15+9),DATA(1:dsf:end,[1]+15+9),...
    DATA(1:dsf:end,[4]+15+9+1),DATA(1:dsf:end,[1]+15+9+1),1e-4);
axis([-150 150 -150 150]), grid on, hold off
axis equal 
legend('estim xy','gps xy');


%% map view 2D (x=E, y=N) with velocity and time and recomputed velocity
figure, plot(DATA(:,[4]+15),DATA(:,[1]+15),'r'), hold on 
plot(DATA(:,[4]+15+9),DATA(:,[1]+15+9),'b') 
% recompute velocity from gps
posN = DATA(:,[1]+15+9);
posE = DATA(:,[4]+15+9);
velN = DATA(:,[1]+15+1+9);
velE = DATA(:,[4]+15+1+9);
%reinterp positions to fill holes and smooth 
deltat = 0.01;
Ti = T(1):deltat:T(end);
posNi = interp1(T,posN,Ti);
posEi = interp1(T,posE,Ti);
posNs = smooth(posNi,100,'loess');
posEs = smooth(posEi,100,'loess');
posNsi = interp1(Ti, posNs,T);
posEsi = interp1(Ti, posEs,T);
%compute velocity vector
velNs = [0; diff(posNs)./deltat];
velEs = [0; diff(posEs)./deltat];
velNsi = interp1(Ti, velNs,T);
velEsi = interp1(Ti, velEs,T);
if 0,
    %verify results: POS
    figure;
    subplot(2,1,1), plot(T,posN,Ti,posNi,Ti,posNs,T,posNsi); legend('original GPS','interpolated GPS','smoothed GPS','re-interpolated GPS');axis([0 T(end),-100,100]);
    subplot(2,1,2), plot(T,posE,Ti,posEi,Ti,posEs,T,posEsi); legend('original GPS','interpolated GPS','smoothed GPS','re-interpolated GPS');axis([0 T(end),-100,100]);
    if exist('T_'),
        subplot(2,1,1), hold on, plot(T_,DATA_(:,[1]+15),'k'); hold off 
        subplot(2,1,2), hold on, plot(T_,DATA_(:,[4]+15),'k'); hold off
    end
    %verify results: VEL
    figure;
    subplot(2,1,1), plot(T,velN,Ti,velNs,T,velNsi); legend('original GPS','smoothed VEL from smoothed POS','re-interpolated');axis([0 T(end),-5,5]);
    subplot(2,1,2), plot(T,velE,Ti,velEs,T,velEsi); legend('original GPS','smoothed VEL from smoothed POS','re-interpolated');axis([0 T(end),-5,5]);
    if exist('T_'),
        subplot(2,1,1), hold on, plot(T_,DATA_(:,[1]+15+1),'k'); hold off 
        subplot(2,1,2), hold on, plot(T_,DATA_(:,[4]+15+1),'k'); hold off
    end
end
%plot vel 
dsf = 50; %down sample factor
quiver(DATA(1:dsf:end,[4]+15),DATA(1:dsf:end,[1]+15),...
    DATA(1:dsf:end,[4]+15+1),DATA(1:dsf:end,[1]+15+1),'g');
quiver(DATA(1:dsf:end,[4]+15+9),DATA(1:dsf:end,[1]+15+9),...
    DATA(1:dsf:end,[4]+15+9+1),DATA(1:dsf:end,[1]+15+9+1),1e-4);
%plot recomputed vel from smoothed pos
quiver(DATA(1:dsf:end,[4]+15),DATA(1:dsf:end,[1]+15),...
    velEsi(1:dsf:end,1),velNsi(1:dsf:end,1),'k');
axis([-150 150 -150 150]), grid on, 
%add time text 
for i= 1:dsf:length(T),
    text(DATA(i,[4]+15),DATA(i,[1]+15),num2str(T(i)));
end
if exist('T_'),
    %simulation was rerun -> plot trajectory and velocity vectors
    % find indexes in T corresponding to times in T_
    T_idxs = [];
    Tsub = T(1:dsf:end);
    for i=1:length(Tsub),
        [v,idx] = min(abs(Tsub(i)-T_));
        T_idxs = [T_idxs, idx];
    end;
    plot(DATA_(:,[4]+15),DATA_(:,[1]+15),'r:','linewidth',3);
    quiver(DATA_(T_idxs ,[4]+15),DATA_(T_idxs,[1]+15),...
        sat(DATA_(T_idxs,[4]+15+1),10),sat(DATA_(T_idxs,[1]+15+1),10),0,'r');
end
hold off
axis equal 
legend('estim xy','gps xy');


%% TRPY 
figure, plot(T, PILOT_COMMANDS), legend('T','R','P','Y');

%command response
figure;
rpy_conv_factor = 15/180*pi;
subplot(3,2,1), plot(T, PILOT_COMMANDS(:,2)*rpy_conv_factor*180/pi,'r'),hold on, plot(T, DATA(:,[2])*180/pi,'b'), hold off
grid on; legend('des roll','roll');
subplot(3,2,2), plot(T, PILOT_COMMANDS(:,3)*rpy_conv_factor*180/pi,'r'),hold on, plot(T, DATA(:,[3])*180/pi,'b'), hold off
grid on; legend('des pitch','pitch');
subplot(3,2,3), plot(T, DATA(:,4+15),'r'), grid on, legend('y');
subplot(3,2,4), plot(T, DATA(:,1+15),'r'), grid on, legend('x');
subplot(3,2,5), plot(T, MODE,'r'), grid on, legend('mode');



%MODE
figure, plot(T, MODE,'r'), grid on, hold off
legend('mode');

%VBAT
figure, plot(T, DATA(:,40),'r'), grid on, hold off
legend('v bat');



%% KF tuning - works after having rerun the filter and saved data into T_ and DATA_
% pos
figure, 
plot(T, DATA(:,[1]+15+9),'r:','linewidth',2), hold on, plot(T, DATA(:,[4]+15+9),'g:','linewidth',2), plot(T, DATA(:,[7]+15+9),'b:','linewidth',2)
plot(T, DATA(:,[1]+15),'r','linewidth',2), plot(T, DATA(:,[4]+15),'g','linewidth',2), plot(T, DATA(:,[7]+15),'b','linewidth',2) 
if exist('T_'),
    plot(T_, DATA_(:,[1]+15),'y','linewidth',3), plot(T_, DATA_(:,[4]+15),'m','linewidth',3), plot(T_, DATA_(:,[7]+15),'c','linewidth',3) 
end
%plot(T_, cumsum(cumsum(DATA_(:,[1]+15+2)-mean(DATA_(:,[1]+15+2)))).*0.02^2,'r.','linewidth',2), ...
%    plot(T_, cumsum(cumsum(DATA_(:,[4]+15+2)-mean(DATA_(:,[4]+15+2)))).*0.02^2,'g.','linewidth',2), ...
%    plot(T_, cumsum(cumsum(DATA_(:,[7]+15+2)-mean(DATA_(:,[7]+15+2)))).*0.02^2,'b.','linewidth',2) 
axis([0 max(T) -150 150]), grid on, hold off
legend('gps x','gps y', 'gps z', 'estim x','estim y','estim z', 'new estim x','new estim y','new estim z');

% vel 
figure, 
plot(T, DATA(:,[1]+15+9+1),'r:'), hold on, plot(T, DATA(:,[4]+15+9+1),'g:'), plot(T, DATA(:,[7]+15+9+1),'b:')
plot(T, DATA(:,[1]+15+1),'r','linewidth',2), plot(T, DATA(:,[4]+15+1),'g','linewidth',2), plot(T, DATA(:,[7]+15+1),'b','linewidth',2) 
plot(T,velNsi,'k');
plot(T,velEsi,'k'); 
if exist('T_'),
    plot(T_, DATA_(:,[1]+15+1),'y','linewidth',3), plot(T_, DATA_(:,[4]+15+1),'m','linewidth',3), plot(T_, DATA_(:,[7]+15+1),'c','linewidth',3) 
    % recompute velocity from currently NEW ESTIMATED position
    posN_ = DATA_(:,[1]+15);
    posE_ = DATA_(:,[4]+15);
    %reinterp positions to fill holes and smooth 
    deltat = 0.01;
    Ti = T_(1):deltat:T_(end);
    posNi_ = interp1(T_,posN_,Ti);
    posEi_ = interp1(T_,posE_,Ti);
    posNs_ = smooth(posNi_,100,'loess');
    posEs_ = smooth(posEi_,100,'loess');
    posNsi_ = interp1(Ti, posNs_,T_);
    posEsi_ = interp1(Ti, posEs_,T_);
    %compute velocity vector
    velNs_ = [0; diff(posNs_)./deltat];
    velEs_ = [0; diff(posEs_)./deltat];
    velNsi_ = interp1(Ti, velNs_,T_);
    velEsi_ = interp1(Ti, velEs_,T_);
    %plot derivative of currently estimated pos
    plot(T_,velNsi_,'k:','linewidth',3.0);
    plot(T_,velEsi_,'k:','linewidth',3.0); 
end
axis([0 max(T) -10 10]), grid on, hold off
legend('gps x','gps y', 'gps z', 'estim x','estim y','estim z', ...
    'vel x from smoothed estim', 'vel y from smoothed estim',...
    'new estim x','new estim y','new estim z');

%acc
figure, 
subplot(3,1,1), plot(T, DATA(:,[1]+15+9+2),'r:'), hold on, 
subplot(3,1,2), plot(T, DATA(:,[4]+15+9+2),'g:'), hold on;
subplot(3,1,3), plot(T, DATA(:,[7]+15+9+2),'b:'), hold on;
subplot(3,1,1), plot(T, DATA(:,[1]+15+2),'r','linewidth',2) 
subplot(3,1,2), plot(T, DATA(:,[4]+15+2),'g','linewidth',2)
subplot(3,1,3), plot(T, DATA(:,[7]+15+2),'b','linewidth',2) 
if exist('T_'),
    subplot(3,1,1), plot(T_, DATA_(:,[1]+15+2),'y','linewidth',3), axis([0 max(T) -10 10]), grid on, hold off
    subplot(3,1,2), plot(T_, DATA_(:,[4]+15+2),'m','linewidth',3), axis([0 max(T) -10 10]), grid on, hold off
    subplot(3,1,3), plot(T_, DATA_(:,[7]+15+2),'c','linewidth',3), axis([0 max(T) -10 10]), grid on, hold off
end
%legend('acc x','acc y', 'acc z', 'estim x','estim y','estim z', 'new estim x','new estim y','new estim z');

%% debug sample loss
figure,
plot(T(2:end),diff(T),'b'), hold on 
plot(T, MODE.*0.1,'g','linewidth',4)
hold off



%% debug position control 
mpc = 10*pi/180; %max_pos_control 
figure; plot(T,DEBUG_POS_CONTROL_x); 
hold on; plot(T,min(max(-mpc,DEBUG_POS_CONTROL_x(:,3)+DEBUG_POS_CONTROL_x(:,4)+DEBUG_POS_CONTROL_x(:,5)),mpc),'b:','linewidth',3); 
plot(T,PILOT_COMMANDS(:,[3]).*15/180*pi);
hold off 
legend('pos ref','pos','kp','ki','kd','PID');
figure; plot(T,DEBUG_POS_CONTROL_y); 
hold on; plot(T,min(max(-mpc,DEBUG_POS_CONTROL_y(:,3)+DEBUG_POS_CONTROL_y(:,4)+DEBUG_POS_CONTROL_y(:,5)),mpc),'b:','linewidth',3); 
plot(T,PILOT_COMMANDS(:,[2]).*15/180*pi);
hold off 
legend('pos ref','pos','kp','ki','kd','PID');

%% debug H control 
%DEBUG_WAYP = [h_ref los los_rate_sat h_wp]
figure, 
plot(T, DEBUG_WAYP);
legend('h ref','los','los rate sat','WP en','h wp');

figure,
plot(T, DEBUG_WAYP(:,1),'b');
hold on;
plot(T, DATA(:,15),'r');
plot(T, DEBUG_WAYP(:,5),'g');
plot(T, -DATA(:,[7]+15+9),'m'), 
plot(T, DEBUG_WAYP(:,4),'c','linewidth',3);
plot(T, MODE,'k');
hold off, grid on
legend('h ref','h','h wp','h baro','wp en','mode');


%% debuf pos control
figure,
plot(T,DEBUG_POS_CONTROL_y(:,1:4))
