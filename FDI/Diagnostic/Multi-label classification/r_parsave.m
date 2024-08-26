function r_parsave(savename,chrname,read,region_len,dark_flag)
    
    save([savename,chrname],'read','region_len','dark_flag','-v7.3')
    
end