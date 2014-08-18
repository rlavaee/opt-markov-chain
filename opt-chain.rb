require 'matrix'

class LehmerCode < Array
  def get_id
    fact = 1
    reverse.each_with_index.inject(0) do |pid,(v,i)|
      pid += v*fact
      fact *= i+1
      pid
    end
  end
end

class Perm < Array

  def get_lcode
    if(@lcode.nil?)
      return @lcode = each_with_index.inject(LehmerCode.new) do |plcode,(v1,i)|
        lcode_i = slice(0,i).inject(v1-1) do |plcode_i,v2|
          (v2 < v1)?(plcode_i-1):(plcode_i)
        end
        plcode << lcode_i
      end
    else
      return @lcode
    end
  end

  def perm_id
    get_lcode.get_id
  end

  def Perm.lcode_to_perm(lcode)
    remaining = (1..lcode.length).to_a
    Perm.new(lcode.inject([]) {|pperm_a,lc| pperm_a << remaining.delete_at(lc)})
  end

  def Perm.id_to_perm(id,len)
    fact = (1..len-1).to_a.inject(:*) || 1
    Perm.lcode_to_perm((0..len-1).to_a.reverse.inject([]) do |pperm_a,i|
      pperm_a << id/fact
      id = id%fact
      fact /= i if(i!=0)
      pperm_a
    end)
  end

  def trim(len)
    return Perm.lcode_to_perm(get_lcode.slice(0,len)+Array.new(length-len,0))
  end

end


class OptState < Perm
  def initialize(*args)
    if(args[0].is_a? Numeric)
      super(Perm.id_to_perm(args[0],args[1]))
    else
      super
    end
  end


  def access_item_at(ind)
    accessed = at(ind)
    osd = nil
    first_move = true
    new_state = OptState.new(each_with_index.inject([length,[]]) do |memo,(v,i)|
      moving = memo[0]
      if(v > accessed)
        memo[1] << v-1
      elsif((v > moving) or first_move)
        osd = i if(first_move)
        memo[1] << moving
        memo[0] = v
        first_move = false
      else
        memo[1] << v
      end
      memo
    end[1])
    return [new_state.perm_id,osd]
  end
end

class TrMatrix < Matrix
  def []=(i, j, x)
    @rows[i][j] = x
  end

  def inspect(m,trim_len=m)
    str = String.new
    fact = (1..m-trim_len).inject(:*) || 1
    str += "\t"
    to_a.each_index do |j|
      str += "#{Perm.id_to_perm(j*fact,m).slice(0,trim_len).inspect}\t"
    end
  
    str+= "\n"
    
    to_a.each_with_index do |matv,j|
      str+="#{Perm.id_to_perm(j*fact,m).slice(0,trim_len).inspect}\t"
      str+= matv.join("\t")+"\n"
    end
    str
  end

  def to_wolfram_str(m,trim_len=m)
    str = String.new
    fact = (1..m-trim_len).inject(:*) || 1
    str+= "{"
  
    column_vectors.to_a.each_with_index do |matv,j|
      str+= "," if (j!=0)
      str+= "{"+matv.to_a.join(",")+"}"
    end

    str+="}"

  end

end

class OptChain
  attr_accessor :t_matrix
  attr_accessor :agg_t_matrix

  def initialize(m)
    @m = m
    @nstates = (1..@m).to_a.inject(:*) || 1
    #@t_matrix = TrMatrix.zero(@nstates)
    @agg_st_size = Array.new
    @n_agg_states = Array.new
    (1..@m).each do |len|
      @agg_st_size[len] = (1..@m-len).inject(:*) || 1
      @n_agg_states[len] = @nstates/@agg_st_size[len]
    end
    compute_agg_matrices
  end


  def compute_agg_st_dist_map(len=@m)
    compute_agg_st_dist(len).row_vectors[0].to_a.each_with_index.inject({}) do |pmap,(p,id)|
      pmap[OptState.new(Perm.id_to_perm(id*@agg_st_size[len],@m)).slice(0,len)] = p
      pmap
    end
  end


  def get_sol_matrix(len=@m)
    p = TrMatrix.columns((@agg_t_matrix[len]-TrMatrix.I(@n_agg_states[len])).column_vectors[0..-2] << TrMatrix.column_vector(Array.new(@n_agg_states[len],1)).column(0))
  end

  def compute_agg_st_dist(len=@m)
    p = TrMatrix.columns((@agg_t_matrix[len]-TrMatrix.I(@n_agg_states[len])).column_vectors[0..-2] << TrMatrix.column_vector(Array.new(@n_agg_states[len],1)).column(0))
    y = TrMatrix.row_vector(Array.new(@n_agg_states[len]-1,0) << 1 )
    y*p.inv
  end

  def compute_agg_matrices
    @agg_t_matrix = Array.new(@m)
    (1..@m).each do |len|
      mat = @agg_t_matrix[len] = TrMatrix.zero(@n_agg_states[len])
      pr = Rational(1,(@m+1)*@agg_st_size[len])

      (0..@nstates-1).each do |state_id|
        state = OptState.new(state_id,@m)
        agg_state_id = state_id/@agg_st_size[len]
        (0..@m-1).each do |ind|
          mat[agg_state_id,state.access_item_at(ind)[0]/@agg_st_size[len]] += pr
        end
        mat[agg_state_id,agg_state_id] += pr
      end
    end
  end

  def compute_osd_dist
    compute_agg_st_dist_map.each.inject([Rational(1,@m+1)]+Array.new(@m,0)) do |posd_dist,(state,p)|
      (0..@m-1).each do |ind|
        posd_dist[state.access_item_at(ind)[1]+1]+=p*Rational(1,@m+1)
      end
      posd_dist
    end
  end

  def compute_osd_dist_v2
    compute_agg_st_dist_map.each.inject([Rational(1,@m+1)]+Array.new(@m,0)) do |posd_dist,(state,p)|
      posd_dist[state.index(@m)+1]+=p*Rational(@m,@m+1)
      posd_dist
    end
  end
end


