require_relative 'opt-chain'
m = ARGV[0].to_i
trim_len = ARGV[1].to_i
chain = OptChain.new(m)
#File.open("osd_dist_#{m}.csv","w") do |f| 
#  f.puts chain.compute_osd_dist.join("\t")
#end

chain.agg_t_matrix.each_with_index do |mat,len|
  if(!mat.nil? and len <= trim_len)
    File.open("chain_#{m}_trim_#{len}.csv","w") do |f|
      f.puts mat.inspect(m,len)
    end
    File.open("chain_wolfram_#{m}_trim_#{len}.csv","w") do |f|
      f.puts chain.get_sol_matrix(len).to_wolfram_str(m,len)
    end
  end
end

(1..trim_len).each do |len|
  File.open("sol_#{m}_trim_#{len}.csv","w") do |f|
    chain.compute_agg_st_dist_map(len).to_a.each do |state,p|
      f.puts "#{state.inspect}\t#{p}\t#{p.to_f}"
    end
  end
end
