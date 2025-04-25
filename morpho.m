function morpho
    %% === Load & Prep ===
    D           = load('P6.mat');
    Fs          = D.sfx;
    dat         = double(D.d);               % [nCh × nTime]
    [nCh,nTime] = size(dat);

    %% === Electrode geometry in μm (1 unit = 10 μm) ===
    nx        = 8;  ny = 4;
    x_elec_um = (0:(nx-1)) * 2 * 50;         % 0→1400 μm
    y_elec_um = (0:(ny-1)) * 2 * 50;         % 0→ 600 μm
    idx_col   = mod(0:nCh-1, nx) + 1;
    idx_row   = floor((0:nCh-1)/nx) + 1;
    x_flat    = x_elec_um(idx_col)';         % 32×1
    y_flat    = y_elec_um(idx_row)';         % 32×1

    %% === 100 Hz analytic signal & phase ===
    [b, afilt] = butter(4, [12 15]/(Fs/2));
    zf          = filtfilt(b, afilt, dat')'; % zero‐phase bandpass
    analytic    = hilbert(zf')';             % complex analytic
    real100     = real(analytic);            % real part [32×nTime]
    phi_wrapped = angle(analytic);           % wrapped phase (–π…π)
    phi100      = unwrap(phi_wrapped, [], 2);% continuous phase along time

    %% === Fine spatial grid ===
    xi_um     = linspace(x_elec_um(1), x_elec_um(end), 700);
    yi_um     = linspace(y_elec_um(1), y_elec_um(end), 300);
    [XI, YI]  = meshgrid(xi_um, yi_um);      % 300×700 grid

    %% === Figure & Axes ===
    hFig = figure('Name','Morpho','Color','w', ...
                  'KeyPressFcn', @keyPress);
    hAx  = axes('Parent', hFig);
    hold(hAx,'on');
    axis(hAx,'vis3d');
    xlim(hAx, [x_elec_um(1) x_elec_um(end)]);
    ylim(hAx, [y_elec_um(1) y_elec_um(end)]);
    zlim(hAx, [-300 300]);
    view(hAx, 45, 30);
    rotate3d(hAx,'on');
    hAx.DataAspectRatioMode    = 'auto';
    hAx.PlotBoxAspectRatioMode = 'manual';
    hAx.PlotBoxAspectRatio     = [7 3 3];

    %% === Pre-create graphics objects ===
    surfObj    = surf(hAx, XI, YI, zeros(size(XI)), ...
                      'EdgeColor','none','FaceAlpha',0.8);
    scatterObj = scatter3(hAx, x_flat, y_flat, real100(:,1), ...
                          75, real100(:,1), 'filled');
    % electrode labels
    loX = 10; loY = 0; loZ = 5;
    lbl = gobjects(nCh,1);
    for i = 1:nCh
        lbl(i) = text(hAx, ...
            x_flat(i)+loX, y_flat(i)+loY, real100(i,1)+loZ, ...
            num2str(i), 'FontSize',8,'FontWeight','bold', ...
            'HorizontalAlignment','left','VerticalAlignment','middle');
    end

    %% === UI Controls ===
    sld = uicontrol('Style','slider', ...
        'Units','normalized','Position',[0.15 0.02 0.7 0.04], ...
        'Min',1,'Max',nTime,'Value',1, ...
        'SliderStep',[1/(2*(nTime-1)), 10/(2*(nTime-1))], ...
        'Callback',@(s,~) updateFrame(s.Value));
    txt = uicontrol('Style','text', ...
        'Units','normalized','Position',[0.88 0.02 0.1 0.04], ...
        'BackgroundColor','w','String','t = 0.000 s');

    %% === Initial Draw ===
    currFrame = 1;
    updateFrame(currFrame);

    %% === Update Function (allows fractional frames) ===
    function updateFrame(fr)
        currFrame = fr;
        t0        = (currFrame-1)/Fs;
        dt        = 1/Fs;

        % fractional interpolation between integer frames
        k0   = floor(currFrame);
        frac = currFrame - k0;
        k1   = min(k0+1, nTime);
        re_t = (1-frac)*real100(:,k0) + frac*real100(:,k1);
        ph_t = (1-frac)*phi100(:,k0)   + frac*phi100(:,k1);

        % temporal slope
        if k0 > 1
            prev_ph = (1-frac)*phi100(:,k0-1) + frac*phi100(:,k0);
            dph_dt  = (ph_t - prev_ph) * Fs;
        else
            dph_dt  = zeros(nCh,1);
        end

        % spatial gradients of phase
        Pmat        = reshape(ph_t, ny, nx);
        [dPdy,dPdx] = gradient(Pmat, 2, 2);
        gX          = dPdx(:);
        gY          = dPdy(:);
        denom       = gX.^2 + gY.^2 + eps;

        % advective velocities
        vx = -dph_dt .* gX ./ denom;
        vy = -dph_dt .* gY ./ denom;

        % advected electrode positions
        x_adv = x_flat + vx*dt;
        y_adv = y_flat + vy*dt;

        % 1) griddata thin‐plate interpolation
        Zraw = griddata(x_adv, y_adv, re_t, XI, YI, 'v4');
        Zraw(isnan(Zraw)) = 0;

        % 2) smooth with 11×11 moving average
        kernel = ones(11,11)/121;
        Zq     = conv2(Zraw, kernel, 'same');

        % update surface
        set(surfObj, 'ZData', Zq, 'CData', Zq);

        % update scatter
        set(scatterObj, ...
            'XData', x_adv, 'YData', y_adv, ...
            'ZData', re_t,  'CData', re_t);

        % update labels
        for i = 1:nCh
            set(lbl(i),'Position',[x_adv(i)+loX, y_adv(i)+loY, re_t(i)+loZ]);
        end

        % update UI
        set(sld, 'Value', currFrame);
        txt.String = sprintf('t = %.3f s', t0);
        title(hAx, sprintf('Real-part 100 Hz, t = %.3f s', t0));

        drawnow;
    end

    %% === Arrow‐key handler ===
    function keyPress(~, evt)
        switch evt.Key
            case 'rightarrow'
                newF = min(currFrame + 0.5, nTime);
            case 'leftarrow'
                newF = max(currFrame - 0.5, 1);
            otherwise
                return;
        end
        updateFrame(newF);
    end
end
