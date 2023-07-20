% gussv907 

uppgift_1a;

uppgift_1b;

function uppgift_1a
    %Uppgift 1A
    
    disp("Uppgift 5.1.1 LAX 1A");
    
    RX = 9369; R2 = 6977; RY = 6294; 
    R4 = 5906; R5 = 1339; R6 = 6243; R7 = 2515;
    E = -14;

    % Matrisen G fås genom att ta "Matrisformen" som har tagits fram i
    % varje nod under nodanalysen 
    G = [ ... 
        -(1/R5+1/RY+1/RX) 1/R5 0; ...           % från nod V1
        1/R5 -(1/R5+1/R6+1/R7) 1/R7; ...        % från nod V2
        0 1/R7 -(1/R2 + 1/R4 + 1/R7) ];         % från nod V3

    J = [ - E/R5 ; E/R5 ; 0 ];

    % V är en vektor som innehåller alla värden för noderna (V1, V2 & V3)
    V = G \ J;      % V = G^-1 * J 

    % Hur vi får fram IX kan tas från uträkningar vid nod V1 (se papper)
    IX = V(1) / RX;
    disp("IX =")
    disp(IX);
end

function uppgift_1b
    %Uppgift 1B
    
    disp("Uppgift 5.1.2 LAX 1B");
    
    RX = 9369; R2 = 6977; RY = 6294; 
    R4 = 5906; R5 = 1339; R6 = 6243; R7 = 2515;
    E = -14;

    J = [ - E/R5 ; E/R5 ; 0 ];      % Samma som tidigare

    RYs = logspace(-3, 10, 10000); % Variera RYs

    IXs = []; % Resultatvektorn

    for ry = RYs

        G = [ -(1/R5+1/ry+1/RX) 1/R5 0; ...           % från nod V1
        
        1/R5 -(1/R5+1/R6+1/R7) 1/R7; ...        % från nod V2
        
        0 1/R7 -(1/R2 + 1/R4 + 1/R7) ];         % från nod V3

        V = G\J;

        IXs(end+1) = V(1) / RX;

    end
    
    figure(1); semilogx(RYs, IXs);
    
    % Extremvärden
    disp("max ström")
    disp(max(IXs))
    disp("min ström")
    disp(min(IXs))
    
    
    PXs = IXs .* IXs * RY; % Beräkna effekten
    figure(2); semilogx(RYs, PXs);
    
    % Hitta var maxeffekten händer
    disp("maxeffekt vid")
    disp(RYs(find(PXs == max(PXs)))) 
end