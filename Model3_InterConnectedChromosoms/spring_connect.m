%% function spring_connect
function [flag_c, dist]=spring_connect(CM, YM, flag_cp, N, dt, l0, koff_0, kon_0)
    dist=sqrt((CM-CM').^2+(YM-YM').^2);
    kon=kon_0*exp(-(dist./l0).^2); %kon decreases with increasing node distances
    koff=koff_0*exp(dist./l0.^2); % koff increases with increasing node distances
    m_on = double(rand(N) < kon*dt); 
    m_off = double(rand(N)< koff*dt & m_on == 0); % m_on and m_off cannot both equal to 1
    flag_c = flag_cp .* (1 - m_off) + (1 - flag_cp) .* m_on;   
    %m_off=0,m_on=1,flg_cp=1,flga_c=1  maintain for establish connection: connection maintained
    %m_off=0,m_on=1,flg_cp=0,flga_c=1  maintain for establish connection: disconnection to connection 
    %m_off=1,m_on=0,flg_cp=1,flga_c=0  disconnect existing connection: connection to disconnection 
    %m_off=1,m_on=0,flg_cp=0,flga_c=0  disconnect existing connection:: disconnecti0on matintained 
    %m_off=0,m_on=0,flg_cp=1,flga_c=1  No change: connection maintained
    %m_off=0,m_on=0,flg_cp=0,flga_c=0  No change:  disconnection maintained
end