
function [ ] = PlotExpData( subplotname, filter_Exp, filter_Ref, color, fluidname )
%PlotExpData Summary of this function goes here
%   Detailed explanation goes here


global Exp Ref PipeVis CFD Rheo

switch subplotname
    case 'p_dpdl_mdot'
        scatter(Exp.MassFlowRate(filter_Exp), Exp.dpdl_p(filter_Exp),'MarkerEdgeColor',color);
        PlotRefData(Ref.MassFlowRate(filter_Ref), Ref.dpdl_p(filter_Ref),color);
        
    case 'a_dpdl_mdot'
        scatter(Exp.MassFlowRate(filter_Exp), Exp.dpdl_ae(filter_Exp),'MarkerEdgeColor',color);
        PlotRefData(Ref.MassFlowRate(filter_Ref), Ref.dpdl_a(filter_Ref),color);
                
    case 'p_f_Usl'
        scatter(Exp.Usl_p(filter_Exp), Exp.f_p(filter_Exp),'MarkerEdgeColor',color);
        PlotRefData(Ref.Usl_p(filter_Ref), Ref.f_p(filter_Ref),color);
                
    case 'a_f_Usl'
        scatter(Exp.Usl_ae(filter_Exp), Exp.f_ae(filter_Exp),'MarkerEdgeColor',color);
        PlotRefData(Ref.Usl_a(filter_Ref), Ref.f_a(filter_Ref),color);
                
    case {'p_f_ReG'}
        scatter(Exp.ReG_p(filter_Exp), Exp.f_p(filter_Exp),'MarkerEdgeColor',color);
        PlotRefData(Ref.ReG_p(filter_Ref), Ref.f_p(filter_Ref),'black');
           
    case {'a_f_ReG'}
        scatter(Exp.ReG_a(filter_Exp), Exp.f_ae(filter_Exp),'MarkerEdgeColor',color);
        PlotRefData(Ref.ReG_a(filter_Ref),Ref.f_a(filter_Ref),'black');
        
    case {'PipeVisco'}
        scatter(PipeVis.NewtSR(filter_Exp), PipeVis.WallSS(filter_Exp),'MarkerEdgeColor',color);
        PlotPipeVisFit(PipeVis.KPrime(filter_Exp), PipeVis.nPrime(filter_Exp), color);
                
end

end
