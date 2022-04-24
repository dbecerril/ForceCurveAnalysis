module ForceCurveAnalysis


using MarketTechnicals
using DelimitedFiles
using Statistics
using Plots
using NoiseRobustDifferentiation
using LsqFit

export SystemParameters, analyzeForceCurveList,analyzeForceCurve,forceCurveOptimizer,writeToFiles
export optimizerMenu
mutable struct SystemParameters
    
    # [FolderPath,MapName,NRows,Ref]
    ReadFilesParam    ::Array
    ElasticConstant   ::Float64 #N/m
    DSensitivity      ::Float64 #V/nm
    ZSensitivity      ::Float64 #nm/step
    OutlierTol        ::Float64
    PointsInWindow    ::Int
    
    # [percentageInit,percentageFinl]
    ZeroPointSubParam ::Array
    
    # [iterations,weigth,dx]
    # default: 200, 0.2, 1
    DervParam         ::Array
    
    #part of the scan you wish to analyze
    #[xi,xf,yi,yf]
    AnalyzeSection    ::Array
    
    # Part of indentation curve usded
    # to calc. young modulus
    ICp                ::Array{Int}
end

function makeFileString(MapName::String, NRows::Int, Ref::String, i::Int,j::Int)
    
    out = string(MapName,"_",i,"_",j,"_",Ref,".txt")
    
    return out

end


function processData_ac(dataIn::Array,SP::SystemParameters)
    
    sigma = SP.OutlierTol        
    nWind = SP.PointsInWindow   
    nps   = SP.ZeroPointSubParam 

    dataIn       = hampelFilter(dataIn, nWind, sigma );
    dataIn       = smoothingAvg(dataIn,nWind);
    dataIn       = zeroLineSubtraction(dataIn,nps[1], nps[2])
    
    return dataIn
    
end


function readFC(path1::String,file_name::String,SP::SystemParameters)
     
    Kc  = SP.ElasticConstant
    DS  = SP.DSensitivity
    ZS  = SP.ZSensitivity
    
    file = string(path1,"/",file_name)
    io   = open(file)
    data = readdlm(io)
    close(io)

    data[:,1] = data[:,1]
    data[:,2] = data[:,2]
    
    return data
    
end

function processData(dataIn::Array,SP::SystemParameters)
    
    sigma = SP.OutlierTol        
    nWind = SP.PointsInWindow   
    nps   = SP.ZeroPointSubParam 
    Kc  = SP.ElasticConstant
    DS  = SP.DSensitivity
    ZS  = SP.ZSensitivity

    dataOut       = zeros(length(dataIn[:,1]),2) 
    dataOut[:,1]  = dataIn[:,1]*ZS
    dataOut[:,2]  = dataIn[:,2]*Kc/DS

    #dataIn       = identifyApproachCurve(dataIn);
    dataOut       = hampelFilter(identifyApproachCurve(dataIn), nWind, sigma );
    dataOut       = smoothingAvg(dataOut,nWind);
    dataOut       = zeroLineSubtraction(dataOut,nps[1], nps[2])
    
    return dataOut
    
end



function hampelFilter(data0::Array, avgWindowSize ::Int64, sigmaTol::Float64 )
    
    nn   = length(data0[:,1])
    di   = Int( round(avgWindowSize/2) )
    outx = []
    outy = []
    
    for i = di + 1 : nn - di
        i1 = i - di
        i2 = i + di   
        
        avgi   = mean(data0[i1:i2,2])
        sigmai = abs( data0[i,2] - avgi ) 
        
        if sigmai < sigmaTol 
            append!(outy,data0[i,2])
            append!(outx,data0[i,1])
        end
        
        
    end
    
    return [outx outy]
end



function smoothingAvg(data0::Array, avgWindowSize ::Int64 )
    
    nn   = length(data0[:,1])
    di   = Int( round(avgWindowSize/2) )
    outx = []
    outy = []
    
    for i = di + 1 : nn - di
        i1 = i - di
        i2 = i + di   
        
        avgi = mean(data0[i1:i2,2])
        append!(outy,avgi)
        append!(outx,data0[i,1])
    end
    
    return [outx outy]
end


####################################################
## pi and pf should be from 1 to 100 with pf > pi
####################################################

function extractCurvePart(data0::Array, pi::Real,pf::Real)
    
    pi > pf ? error("parameter 2 (pinitial) must be larger than parameter 1 (pfinal)") : pi = pi
    
    nn = length(data0[:,1])
    ni = Int(round(nn*pi/100) )
    
    ni == 0 ? ni = 1 : ni = ni
    
    nf = Int(round(nn*pf/100) )
    
    dataOut = data0[ni:nf,:]

    return dataOut
end







function zeroLineSubtraction(data0::Array{Float64}, fcInterval0::Real,fcIntervalf::Real)
    
    data1 = extractCurvePart(data0::Array, fcInterval0 ,fcIntervalf)
    
    model(x,p) = p[1]*x .+ p[2] 
    p0         = [0.5,0.5]

    fit = curve_fit(model, data1[:,1], data1[:,2], p0)
    pp = fit.param

    return data0 .- pp[2]
    
    
end

function approachCurveLinearFit(ac::Array{Float64}, acInterval0::Real,acIntervalf::Real)
    data1 = extractCurvePart(ac::Array, acInterval0 ,acIntervalf)
    model(x,p) = p[1]*x .+ p[2] 
    p0         = [0.1,0.1]

    fit = curve_fit(model, data1[:,1], data1[:,2], p0)
    pp = fit.param
    
    return  pp
    
    
end





function identifyApproachCurve(data0::Array)

    Indx = findmax(data0[:,1])[2]
    datatemp = data0[1:Indx,:]
    Indx2 = findmax(datatemp[:,2])[2]

    datatemp = datatemp[1:Indx2,:]

    return datatemp
end

function identifyApproachCurve(data0::Array,cpIndx::Int)
    
    NN = length(data0[:,1])

    return data0[cpIndx:NN,:]
end


function findFirstDerivativeZero(dervData0::Array)
    
    NN      = length(dervData0)
    indxOut = 1
    
    for i  = 1:NN-1
        i1 = NN - i + 1
        i2 = NN - i
        
        if dervData0[i2] <= 0 <= dervData0[i1]
            indxOut = i1
            break
        end
        
    end
    
    return indxOut
end


function findContactPointIndx(data0::Array,Niter::Real,alpha::Real,dx1::Real)

    
    #u     = TVRegDiff(data0[:,2], Int(Niter), alpha, dx = dx1)
    u     = tvdiff(data0[:,2], Int(Niter), alpha, dx = dx1)

    indx0 = findFirstDerivativeZero(u)    
    
    return indx0
end


function makeIndentationCurve(acData0::Array{Float64}, SP::SystemParameters)
    
    Kc = SP.ElasticConstant
    S  = SP.DSensitivity
    acData0[:,1] = acData0[:,1] .- acData0[1,1]
    m = Kc
    
    outx = []
    outy = []
    
    for i = 1:length(acData0[:,1])
        deltai =  ( acData0[i,1] - (acData0[i,2]-acData0[1,2])/m   )^2
        #deltai =  ( acData0[i,1] - (acData0[i,2])/m   )^2

        fi     = acData0[i,2]
        
        outx = append!(outx,deltai)
        outy = append!(outy,fi)
        
    end
      
    return [outx outy]
end

function indentationCurveLinearFit(ac::Array{Float64}, acInterval0::Real,acIntervalf::Real)
    data1 = extractCurvePart(ac::Array, acInterval0 ,acIntervalf)
    
    model(x,p) = p[1]*x .+ p[2] 
    p0         = [minimum(data1[:,2]),minimum(data1[:,2])]

    fit = curve_fit(model, data1[:,1], data1[:,2], p0)
    pp = fit.param
    
    return  pp
    
    
end

###
## acData0: Must be only the approach curve, not the full curve.
###
function youngModulusFromIC(acData0::Array,SP::SystemParameters)
    
    pi = SP.ICp[1]
    pf = SP.ICp[2]
    nu = 0.5
    C0 = 1
    
    ic  = makeIndentationCurve(acData0, SP)
    
    if length(ic[:,1]) > 4 
        pps = indentationCurveLinearFit(ic, pi, pf)
    else
        pps= [0,0]
    end
    
    tanalpha = 0.4348123749609336 # tan(alpha)
    term1 =   ( 3*pi*tanalpha  /(8*C0) )

    E = term1*pps[1]
    
    return (E,pps)
end



function makeIndentationCurveFromParabola(acData0::Array{Float64}, SP::SystemParameters)
    
    Kc = SP.ElasticConstant
    S  = SP.DSensitivity

    acData0[:,1] = acData0[:,1] .- acData0[1,1]
    yh     = Kc*acData0[:,1]    
    
    Kc = SP.ElasticConstant
    pp           = approachCurveParabolaFit(acData0)
    dd = (yh/Kc - ( real.( sqrt.(Complex.(- 4*pp[1]*yh.+pp[2]^2 ) ) ) .+ pp[2]  )/(2*pp[1] ) ).^2
    return [dd yh]
end

function approachCurveParabolaFit(acData )
    
    model(x,p) = p[1]*x.^2 .+ p[2]*x

    p0         = [minimum(acData[:,2]),minimum(acData[:,2])]

    fit = curve_fit(model, acData[:,1], acData[:,2], p0)
    pp = fit.param
    
    return  pp
    
    
end


function younModulusFromParabolaFit(acData,SP::SystemParameters)
   
    ic = makeIndentationCurveFromParabola(acData::Array{Float64}, SP::SystemParameters)
    pi = SP.ICp[1]
    pf = SP.ICp[2]
    alpha = 23.5*(pi/180)
    C0 = 1
    coef = 3*pi*tan(alpha)/(8*C0)

    if length(ic[:,1]) > 4 
        pps = indentationCurveLinearFit(ic, pi, pf)
    else
        pps= [0,0]
    end

    E = ( 3*pi*tan(alpha)/(8*C0) )*pps[1]
    
    
    return E

end


function analyzeForceCurve(SP::SystemParameters)
    
    mapNxNy   = SP.AnalyzeSection
    NameParam = SP.ReadFilesParam
    Dparam    =  SP.DervParam
    counter   = 1
   
    mapNxNy[1] == 0 ? mapNxNy =[1,NameParam[3],1,NameParam[3] ]  : mapNxNy = mapNxNy

    Dims      = ( mapNxNy[2]-mapNxNy[1] + 1 ) * ( mapNxNy[4]-mapNxNy[3] + 1 )
    out       = zeros(Dims,5)

    for  j = mapNxNy[3]:mapNxNy[4], i = mapNxNy[1]:mapNxNy[2]
        
    
        file_ij = makeFileString(NameParam[2], NameParam[3], NameParam[4], i,j)
        data_ij = readFC(NameParam[1],file_ij,SP)
        data_ij = processData(data_ij,SP)
        cpIndx  = findContactPointIndx(data_ij,Dparam[1],Dparam[2],Dparam[3])
        acData  = identifyApproachCurve(data_ij,cpIndx)
                
        
        #young_poli = younModulusFromParabolaFit(acData,SP::SystemParameters)
        ppi     = approachCurveLinearFit(acData, SP.ICp[1],SP.ICp[2])
        slopei  = ppi[1]
        
        # Some expection handling
        # This occurs in "calcYoungModulus" when the ic has too few points
        
        youngi  = youngModulusFromIC(acData::Array,SP::SystemParameters)[1]     
        youngi == 0 ? println("ic $(file_ij) has too few points youngMod = 0") : j = j
      
        
        out[counter,1 ] = i
        out[counter,2 ] = j

        
        #contac point
        out[counter,3 ] = 65000 - data_ij[cpIndx,1]
        out[counter,4 ] = slopei#young_poli 
        out[counter,5]  = youngi
        
        counter +=1
        println("point $i $j complete")
    end
    
    return out

end







function forceCurveOptimizer(ListPoints::Array{Int64} ,show ::Int64 , SP::SystemParameters)
    
    mapNxNy   = SP.AnalyzeSection
    NameParam = SP.ReadFilesParam
    Dparam    =  SP.DervParam
    Dims      = ( mapNxNy[2]-mapNxNy[1] + 1 ) * ( mapNxNy[4]-mapNxNy[3] + 1 )
    Nlist     = length(ListPoints[:,1])
    global out       = []
    
    for  i = 1:Nlist
        
        xi = ListPoints[i,1]
        yi = ListPoints[i,2]
        
        file_ij  = makeFileString(NameParam[2], NameParam[3], NameParam[4], xi,yi)
        data_ij0 = readFC(NameParam[1],file_ij,SP)
        data_ij0 = hampelFilter(data_ij0, 4, SP.OutlierTol )

        data_ij = processData(data_ij0,SP)
        
        cpIndx  = findContactPointIndx(data_ij,Dparam[1],Dparam[2],Dparam[3])
        
        acData  = identifyApproachCurve(data_ij,cpIndx)
        
        
        if show == 1
           aploti  = plot(data_ij0[:,1],data_ij0[:,2],title = " x = $(xi), y = $(yi)", marker = :dot,
                        frame_style = :box, label = "force curve")
                      #vline!([data_ij[cpIndx,1]],label = "contact point") 
                      #xlims!(data_ij[Int(round( (length(data_ij[:,1])/2 ))),1],data_ij[length(data_ij[:,1]),1])
             out = push!(out,aploti)
            @show maximum(data_ij0[:,2])
        elseif show == 2
            aploti  = plot(data_ij[:,1],data_ij[:,2],title = " x = $(xi), y = $(yi)", marker = :dot,
                        frame_style = :box, label = "force curve")
                vline!([data_ij[cpIndx,1]],label = "contact point") 
             #xlims!(data_ij[Int(round( (length(data_ij[:,1])/2 ))),1],data_ij[length(data_ij[:,1]),1])
             out = push!(out,aploti)
        elseif show == 3
            #pps     = approachCurveParabolaFit(acData )
            ppi     = approachCurveLinearFit(acData, 1,100)

            @show slopei  = ppi[1]
            bi      = ppi[2]

            aploti  = plot(acData[:,1].-acData[1,1],acData[:,2],title = " x = $(xi), y = $(yi) ", marker = :dot,
                      frame_style = :box, label = "ac")
                    
                    plot!(acData[:,1].-acData[1,1],slopei*(acData[:,1].-acData[1,1]) , label = "linear fit", lw = 2)
                    #plot!(acData[:,1].-acData[1,1],pps[1]*acData[:,1].^2 + pps[2]*acData[:,1], label = "parabolic fit", lw = 1 )           
                    plot!(acData[:,1].-acData[1,1],SP.ElasticConstant*(acData[:,1].-acData[1,1]), label = "hard surf.", lw = 2 )           
                    out = push!(out,aploti)

    #    elseif show == 4
    #        out3 = younModulusFromParabolaFit(acData,SP::SystemParameters)

    #        ic   = makeIndentationCurveFromParabola(acData::Array{Float64}, SP::SystemParameters)

    #        aploti  = plot(ic[:,1],ic[:,2],title = " x = $(xi), y = $(yi)", marker = :dot,
    #                  frame_style = :box, label = "ic from parab. fit" )
                   
    #                 out = push!(out,aploti)
 
        elseif show == 4
            ic  = makeIndentationCurve(acData, SP)
            #ic      = smoothingAvg(ic,5)
            nn = length(ic[:,1])
            ni = Int(round(nn*SP.ICp[1]/100) )
            nf = Int(round(nn*SP.ICp[2]/100) ) 
            
            youngi  = youngModulusFromIC(acData,SP) 
            mm      = youngi[2]    
            @show youngi[1]
            aploti  = plot(ic[:,1],ic[:,2],title = " x = $(xi), y = $(yi) ", marker = :dot,
                        frame_style = :box, label = "Indentation curve")
                    plot!(ic[:,1],mm[1]*ic[:,1] .+ mm[2],label = "Section l.fit")
            vline!([ic[ni+1,1],ic[nf,1]],label = :none)
            out = push!(out,aploti)

        end
        
        
    end
    
    return out
end

function writeToFiles(data::Array,date::String,SP::SystemParameters)
    rfp = SP.ReadFilesParam 
    icp = SP.ICp

    output_path  =  string(rfp[1],"/") #string(pwd({}),"/outputfiles/")
    output_name  = string(output_path, rfp[2],"_",date,"_ic$(icp[1])_$(icp[2])_topography.xyz")
    output_name2  = string(output_path,rfp[2],"_",date,"_ic$(icp[1])_$(icp[2])_ac_linearfit.xyz")
    output_name3  = string(output_path,rfp[2],"_",date,"_ic$(icp[1])_$(icp[2])_young_ic_linearfit.xyz")
    
    io = open(output_name,"w")
    writedlm(io,[data[:,1:2] data[:,3]])
    close(io)
    
    io = open(output_name2,"w")
    writedlm(io,[data[:,1:2] data[:,4]])
    close(io)
    
    io = open(output_name3,"w")
    writedlm(io,[data[:,1:2] data[:,5] ] )
    close(io)
end


function plotListPoint(plotArray)
    
    nplots = length(plotArray)
    
    nn  = Int( round(nplots/2) )
    
    ap = plot(collect(plotArray), layout = (nn,nn),size = (500,500) )
    
    display(ap)
    
end





function optimizerMenu()

    println("Options for the show variable are:")
    println("show = 1 : full approach curve")
    println("show = 2 : AC with contact point. Zoomed near CP")
    println("show = 3 : AC with parabolic and Linear fit")
    println("show = 4 : Indentation curve from parabolic fit")
    println("show = 5 : Indentation curve made directly from AC")

end


function analyzeForceCurveList(FileNameList::Array{String},path0::String ,show ::Int64 , SP::SystemParameters)
    
    mapNxNy   = SP.AnalyzeSection
    NameParam = SP.ReadFilesParam
    Dparam    =  SP.DervParam
    Dims      = ( mapNxNy[2]-mapNxNy[1] + 1 ) * ( mapNxNy[4]-mapNxNy[3] + 1 )
    Nlist     = length(FileNameList[:,1])
    out       = []
    
    for  i = 1:Nlist
        
        file_i  = FileNameList[i] 
        
        data_ij = readFC(path0,file_i,SP)
        data_ij[:,2] = reverse(data_ij[:,2])
       # data_ij[:,1] =  -(data_ij[:,1] .- 6000) 
        data_ij = processData_ac(data_ij,SP)
        
        cpIndx  = findContactPointIndx(data_ij,Dparam[1],Dparam[2],Dparam[3])
        
        acData  = identifyApproachCurve(data_ij,cpIndx)

        ppi     = approachCurveLinearFit(acData, 1,100)

        slopei  = ppi[1]
        bi     = ppi[2]


        #

      
        if show == 1
            aploti  = plot(data_ij[:,1],data_ij[:,2],title = "$(FileNameList[i])", marker = :dot,
                        frame_style = :box, label = "force curve")
                  vline!([data_ij[cpIndx,1]],label = "contact point") 
                  #xlims!(data_ij[Int(round( (length(data_ij[:,1])/2 ))),1],data_ij[length(data_ij[:,1]),1])
            
        elseif show == 2
            aploti  = plot(acData[:,1],acData[:,2],title = "$(FileNameList[i]) : \n ac compared with fit of curve", marker = :dot,
                      frame_style = :box, label = "aprch. curve")
                      plot!(acData[:,1],slopei*acData[:,1].+bi, label = "m = $(slopei)")
            
        elseif show == 3
            aploti  = plot(acData[:,1].-acData[1,1],acData[:,2],title = "$(FileNameList[i])", marker = :dot,
                      frame_style = :box, label = "ac compared to hard surface ac")
                      plot!(acData[:,1].-acData[1,1],SP.ElasticConstant*(acData[:,1].-acData[1,1]) )           
            
        elseif show == 4
            ic      = makeIndentationCurve(acData, SP)
            ic      = smoothingAvg(ic,5)
            nn = length(ic[:,1])
            ni = Int(round(nn*SP.ICp[1]/100) )
            nf = Int(round(nn*SP.ICp[2]/100) ) 
            
            youngi  = calcYoungModulus(acData,SP)     
            @show youngi
            aploti  = plot(ic[:,1],ic[:,2],title = "$(FileNameList[i]) ", marker = :dot,
                        frame_style = :box, label = "Indentation curve")
            vline!([ic[ni,1],ic[nf,1]],label = :none)


        end
        
        out = push!(out,aploti)
        
    end
    
    return out
end



end