require_relative 'opt-chain'
require 'test/unit'

class TestPerm < Test::Unit::TestCase
  def test_lehmer_code 
    assert_equal([0,0,0],Perm.new([1,2,3]).to_lehmer_code)
    assert_equal([2,1,0],Perm.new([3,2,1]).to_lehmer_code)
    assert_equal([1,2,1,0],Perm.new([2,4,3,1]).to_lehmer_code)
    assert_equal([1,2,2,0,0],Perm.new([2,4,5,1,3]).to_lehmer_code)
  end

  def test_lehmer_code_rev
    assert_equal([1,2,3],Perm.lcode_to_perm([0,0,0]))
    assert_equal([3,2,1],Perm.lcode_to_perm([2,1,0]))
    assert_equal([2,4,3,1],Perm.lcode_to_perm([1,2,1,0]))
    assert_equal([2,4,5,1,3],Perm.lcode_to_perm([1,2,2,0,0]))
  end

  def test_perm_id 
    assert_equal(0,Perm.new([1,2,3]).perm_id)
    assert_equal(5,Perm.new([3,2,1]).perm_id)
    assert_equal(11,Perm.new([2,4,3,1]).perm_id)
    assert_equal(40,Perm.new([2,4,5,1,3]).perm_id)
  end

  def test_perm_id_rev 
    assert_equal([1,2,3],Perm.id_to_perm(0,3))
    assert_equal([3,2,1],Perm.id_to_perm(5,3))
    assert_equal([2,4,3,1],Perm.id_to_perm(11,4))
    assert_equal([2,4,5,1,3],Perm.id_to_perm(40,5))
  end
end

class TestOptChain < Test::Unit::TestCase
  def test_transition
    assert_equal([5,3,4,1,2],Perm.id_to_perm(OptState.new([2,4,5,1,3]).access_item_at(4),5))
    assert_equal([5,2,4,1,3],Perm.id_to_perm(OptState.new([2,4,5,1,3]).access_item_at(2),5))
    assert_equal([6,1,2,3,5,4],Perm.id_to_perm(OptState.new([1,2,3,4,6,5]).access_item_at(3),6))
    assert_equal([10,8,9,1,2,3,4,5,6,7],Perm.id_to_perm(OptState.new([8,9,10,1,2,3,4,5,6,7]).access_item_at(2),10))

  end

  def test_chain
    (2..6).to_a.each do |m|
      chain = OptChain.new(m)
      dist = chain.compute_st_dist
      puts "\nDistribution is: #{dist}"
      assert_equal(dist,dist*chain.t_matrix)
    end
  end
end

