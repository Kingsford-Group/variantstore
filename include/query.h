/*
 * ============================================================================
 *
 *       Filename:  query.h
 *
 *         Author:  Prashant Pandey <ppandey2@cs.cmu.edu>
 *									Yinjie Gao <yinjieg@andrew.cmu.edu>
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef __QUERY_H__
#define __QUERY_H__

#include <string>
#include <vector>
#include <map>
#include <cassert>
#include <iostream>
#include <fstream>
#include "index.h"
#include "variant_graph.h"

std::string REF = "ref";

namespace variantstore {

	struct Variant {
		uint64_t var_pos;
		std::string ref;
		std::string alt;
		// sample names and phasing info.
		std::vector<std::pair<std::string, std::string>> samples;
	};

	void print_var (Variant *var, ofstream& out) {
		out << var->var_pos << "\t"<< var->ref << "\t" << var->alt << "\t";
		for (auto sample : var->samples) {
			out << sample.first << "(" << sample.second << ") ";
		}
		out << std::endl;
		return;
	}


	/* ----------------------------------------------------------------------------
		 Support Func: traverse back until a node containg sample_id is found
		 return such node, record the pos in ref and sample
		 */
	Graph::vertex get_prev_vertex_with_sample ( VariantGraph *vg, Index *idx,
																							const uint64_t pos,
																							const string sample_id,
																							uint64_t &ref_pos,
																							uint64_t &sample_pos)
	{
		// the node right before cur_pos
		uint64_t cur_pos = pos;
		uint64_t cur_ref_node_idx;
		Graph::vertex v = idx->find(cur_pos, cur_ref_node_idx);
		Graph::vertex v_find = v;
		VariantGraphVertex::sample_info sample;
		VariantGraphVertex::sample_info sample_find;
		bool sample_found = false;

		while (true) { // do until find sample or to the front of graph
			v = idx->previous(cur_ref_node_idx);

			console->debug("cur_ref_node_idx: {}", cur_ref_node_idx);
			// if reach the first node in the graph, return the first node
			if (cur_ref_node_idx <= 1) {
				ref_pos = 1;
				v_find = v;
				assert(v_find == 0);
				vg->get_sample_from_vertex_if_exists(v_find, REF, sample_find);
				sample_pos = sample_find.index();
				break;
			}

			// check all outgoing edges of v by BFS, report of a node containing sample
			// is found. ref_pos is updated whenever a node with ref is encountered
			VariantGraph::VariantGraphIterator it = vg->find(v, 1);
			++it;
			while(!it.done()) {
				v = (*it)->vertex_id();

				// TODO: consider consecutive mutation
				if (vg->get_sample_from_vertex_if_exists(v, REF, sample))
				{ ref_pos = sample.index(); }
				if (vg->get_sample_from_vertex_if_exists(v, sample_id, sample_find))
				{
					v_find = v;
					sample_found = true;
					sample_pos = sample_find.index();
				}
				++it;
				cur_ref_node_idx--;
			}

			if (sample_found == true) {break;}
		}

		console->debug("The closest node before {} contain sample {} is node {}. Ref Index: {}. Sample index: {}",
									 pos, sample_id, v_find, ref_pos, sample_pos);

		return v_find;
	} // get_prev_vertex_with_sample()


	/* ----------------------------------------------------------------------------
		 Give positional range and sample id
		 return sequence of sample_id in ref position [pos_x, pos_y)
		 */
	std::string query_sample_from_ref ( VariantGraph *vg, Index *idx,
																			const uint64_t pos_x, const uint64_t pos_y,
																			const string sample_id,
																			bool print=false, std::string outfile="")
	{
		std::string seq = "";
		uint64_t ref_pos;
		uint64_t sample_pos;
		// get the closest node containing sample_id upstream pos_x
		Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
																													sample_id,
																													ref_pos, sample_pos);

		// Do DFS on the path of the sample record seq from pos_x to pos_y in ref
		VariantGraph::VariantGraphPathIterator it = vg->find(closest_v, sample_id);
		bool record_seq = false;
		std::string temp;

		while (!it.done()) {
			temp.assign(vg->get_sequence(*(*it))); // store seq of current node
			console->debug("Check node: {}, reference Index: {}, seq: {}",
										 (*it)->vertex_id(), ref_pos, temp);
			uint64_t l = (*it)->length();

			// if there no outgoing node containing ref
			// TODO: consider consecutive mutation
			uint64_t next_ref_pos = ref_pos + l;

			// Check outgoing nodes, update ref_pos to the index of node containing ref
			VariantGraph::VariantGraphIterator bfs_it = vg->find((*it)->vertex_id(), 1);
			++bfs_it;
			while (!bfs_it.done()) {
				Graph::vertex v = (*bfs_it)->vertex_id();
				VariantGraphVertex::sample_info sample;
				if (vg->get_sample_from_vertex_if_exists(v, REF, sample)){
					next_ref_pos = sample.index();
					break;
				}
				++bfs_it;
			}

			console->debug("Next reference index: {}", next_ref_pos);

			if (record_seq == true && next_ref_pos < pos_y) {
				// already start recoding and no reach pos_y at this node
				seq += temp;
			} else if (record_seq == true && next_ref_pos >= pos_y) {
				// already start recoding and no reach pos_y at this node
				seq += temp.substr(0, pos_y - ref_pos);
				break;
			} else if (next_ref_pos >= pos_x && next_ref_pos < pos_y) {
				// reach pos_x at this node and not reach pos_y
				record_seq = true;
				seq += temp.substr(pos_x - ref_pos);
			} else if (next_ref_pos >= pos_x && next_ref_pos >= pos_y) {
				// This node covers [pos_x, pos_y)
				seq = temp.substr(pos_x - ref_pos, pos_y - pos_x);
				break;
			}

			++it;
			ref_pos = next_ref_pos; // update ref_pos
		}

		if (print == true)
		{
			ofstream out;
			out.open(outfile);
			out << seq << std::endl;
			out.close();
		}

		return seq;
	} // query_sample_from_ref()

	/* ----------------------------------------------------------------------------
		 Give positional range and sample id
		 return sequence of sample_id in sample's position [pos_x, pos_y)
		 */
	std::string query_sample_from_sample ( VariantGraph *vg, Index *idx,
																				 const uint64_t pos_x,
																				 const uint64_t pos_y,
																				 const string sample_id,
																				 bool print=false, std::string outfile="")
	{
		std::string seq = "";
		uint64_t ref_pos;
		uint64_t sample_pos;

		Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
																													sample_id,
																													ref_pos, sample_pos);
		// look for prev node containing sample until such node's sample pos < pos_x
		while (sample_pos >= pos_x && closest_v > 0)
		{
			uint64_t pos = ref_pos;
			closest_v = get_prev_vertex_with_sample(vg, idx, pos, sample_id,
																							ref_pos, sample_pos);
		}

		console->debug("Closest notes containing sample with sample_pos < {}: {}",
									 pos_x, closest_v);

		// DFS from such node and record seq btw [pos_x, pos_y) in sample's coordinate
		VariantGraph::VariantGraphPathIterator it = vg->find(closest_v, sample_id);
		bool record_seq = false;
		std::string temp;

		while (!it.done()) {
			temp.assign(vg->get_sequence(*(*it)));
			console->debug("Check node: {}, sample Index: {}, seq {}",
										 (*it)->vertex_id(), sample_pos, temp);
			uint64_t l = (*it)->length();
			uint64_t next_sample_pos = sample_pos + l;

			if (record_seq == true && next_sample_pos < pos_y) {
				seq += temp;
			} else if (record_seq == true && next_sample_pos >= pos_y) {
				seq += temp.substr(0, pos_y - sample_pos);
				break;
			} else if (next_sample_pos >= pos_x && next_sample_pos < pos_y) {
				record_seq = true;
				seq += temp.substr(pos_x - sample_pos);
			} else if (next_sample_pos >= pos_x && next_sample_pos >= pos_y) {
				seq = temp.substr(pos_x - sample_pos, pos_y - pos_x);
				break;
			}

			++it;
			sample_pos = next_sample_pos;
		}

		if (print == true)
		{
			ofstream out;
			out.open(outfile);
			out << seq << std::endl;
			out.close();
		}

		return seq;
	} // query_sample_from_sample()


	/* ----------------------------------------------------------------------------
		 Support Func: get samples in a nodes
		 Return true if the node contains samples other than ref
		 */
	bool get_samples(const VariantGraphVertex* v, VariantGraph *vg,
									 std::vector<std::pair<std::string, std::string>>
									 &sample_ids)
	{
		bool is_var = false;
		sample_ids = {};
		for (int i = 0; i < v->s_info_size(); ++i) {
			//const VariantGraphVertex::sample_info& s = v->s_info(i);
			std::string sample_id =
				vg->get_sample_name(vg->get_sample_id(*v, i));
			if (sample_id != REF) {
				std::string phasing = vg->get_sample_phasing(*v, i); 					
				sample_ids.push_back(std::make_pair(sample_id, phasing));
				is_var = true;
			}
		}
		return is_var;
	}





	/* ----------------------------------------------------------------------------
		 Support Func: Return true if found the closest variants at or after the
		 given pos. Return false if no variant is found reaching the end of the
		 graph
		 */

	bool next_variant_in_ref ( VariantGraph *vg, Index *idx,
														 const uint64_t pos,
														 vector<Variant> &vars,
														 uint64_t &next_pos)
	{
		bool found_var = false;
		Graph::vertex v = idx->find(pos); // the node before pos
		VariantGraph::VariantGraphPathIterator it = vg->find(v, REF);
		VariantGraph::VariantGraphPathIterator next_it = vg->find(v, REF);
		++next_it;

		while (!it.done()) {
			// console->debug("Check vertex: {}", (*it)->vertex_id());
			VariantGraph::VariantGraphIterator bfs_it=vg->find((*it)->vertex_id(), 1);
			++bfs_it;
			uint32_t num_out_nodes = 0;
			// check all outgoing edges and record alternatives
			while (!bfs_it.done())
			{
				Variant var;
				if ((*bfs_it)->vertex_id() == (*next_it)->vertex_id()) {
					++bfs_it;
					continue;
				}
				// console->debug("Check bfs vertex: {}", (*bfs_it)->vertex_id());
				num_out_nodes++;
				//Graph::vertex out_v = (*bfs_it)->vertex_id();
				std::vector<std::pair<std::string, std::string>> sample_ids;
				// got variant nodes' samples_ids
				if (get_samples((*bfs_it), vg, sample_ids)) { // found var
					VariantGraphVertex::sample_info sample;

					// Deletion
					if (vg->get_sample_from_vertex_if_exists((*bfs_it)-> vertex_id(),
																									 REF, sample)) {
						console->debug("Found a deletion");
						std::string alt = "";
						std::string ref = vg->get_sequence(*(*next_it));
						var.ref.assign(ref);
						var.alt.assign(alt);
						var.samples = sample_ids;
						if (vg->get_sample_from_vertex_if_exists((*next_it)-> vertex_id(),
																										 REF, sample)) {
							var.var_pos = sample.index();
						} else {
							console->error("{} is not in the reference path", (*next_it)->
														 vertex_id());
						}
					} else {
						vg->get_sample_from_vertex_if_exists((*it)-> vertex_id(), REF,
																								 sample);
						uint64_t prev_ref_idx = sample.index();
						// dfs from bfs_it in sample's path
						VariantGraph::VariantGraphPathIterator dfs_it =
							vg->find((*bfs_it)->vertex_id(), sample_ids[0].first);
						++dfs_it;
						if (!vg->get_sample_from_vertex_if_exists((*dfs_it)-> vertex_id(),
																											REF, sample))
						{
							console->error("consecutive mutation near {}", pos);
						}

						uint64_t next_ref_idx = sample.index();
						std::string prev_ref = vg->get_sequence(*(*it));

						// Insertion
						if (next_ref_idx == prev_ref_idx + prev_ref.length()) {
							// console->debug("Found an insertion");
							std::string ref = "";
							var.ref.assign(ref);
							std::string alt = vg->get_sequence(*(*bfs_it));
							var.alt.assign(alt);
							var.samples = sample_ids;
							var.var_pos = next_ref_idx - 1;
						} else {
							// substitution
							std::string alt = vg->get_sequence(*(*bfs_it));
							var.alt.assign(alt);
							std::string ref = vg->get_sequence(*(*next_it));
							var.ref.assign(ref);
							var.samples = sample_ids;
							if (vg->get_sample_from_vertex_if_exists((*next_it)-> vertex_id(),
																											 REF, sample)) {
								var.var_pos = sample.index();
							} else {
								console->error("{} is not in the reference path", (*next_it)->
															 vertex_id());
							}
						}
					}
				}

				// only add var if not seen before.
				if (vars.size() < 1 || (vars.back().var_pos != var.var_pos)) {
					found_var = true;
					vars.push_back(var);
				}
				++bfs_it;
			}

			if (found_var == true) {
				break;
			}

			++it;
			++next_it;
		} // end DFS in ref

		VariantGraphVertex::sample_info sample;
		if (vg->get_sample_from_vertex_if_exists((*next_it)-> vertex_id(),
																						 REF, sample)) {
			next_pos = sample.index();
		} else {
			console->error("{} is not in the reference path", (*next_it)->
										 vertex_id());
		}

		return found_var;
	}

	/* ----------------------------------------------------------------------------
		 Find a variant near given position
		 */
	bool closest_var ( VariantGraph *vg, Index *idx, const uint64_t pos,
										 std::vector<Variant> &vars, bool print=false, std::string
										 outfile="")
	{
		std::vector<Variant> next_var;
		uint64_t next_pos;

		if (next_variant_in_ref(vg, idx, pos, next_var, next_pos)) {
			uint64_t next_var_pos = next_var[0].var_pos;
			std::vector<Variant> prev_var;
			int cur_pos = pos-(next_var_pos-pos);
			if (cur_pos > 0) {
				next_variant_in_ref(vg, idx, cur_pos, prev_var, next_pos);
				uint64_t prev_var_pos = prev_var[0].var_pos;
				if (prev_var_pos != next_var_pos) {
					vars = prev_var;
				} else {
					vars = next_var;
				}
			} else {
				vars = next_var;
			}
		} else {
			int cur_pos = pos - 1;
			while (cur_pos > 0 && !next_variant_in_ref(vg, idx, cur_pos, next_var,
																								 next_pos)) {
				if (cur_pos == 1) { return false;}
				cur_pos--;
			}
			vars = next_var;
		}

		if (print==true) {
			ofstream out;
			out.open(outfile);
			for (auto var : vars)
				print_var(&var, out);
			out.close();
		}

		return true;
	}


	/* ----------------------------------------------------------------------------
		 return all variants in sample_id occuring in sample's coordinate
		 [pos_x, pos_y)
		 */
	std::vector <Variant> get_sample_var_in_sample ( VariantGraph *vg, Index *idx,
																									 const uint64_t pos_x,
																									 const uint64_t pos_y,
																									 const std::string sample_id,
																									 bool print=false, std::string outfile="")
	{
		// traverse in sample's coordinate from pos_x to pos_y
		// report if the sample's node is a variant node
		std::vector <Variant> vars;
		uint64_t ref_pos;
		uint64_t sample_pos;
		Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
																													sample_id,
																													ref_pos, sample_pos);

		// look for prev node containing sample until such node's sample pos < pos_x
		while (sample_pos >= pos_x && closest_v > 0)
		{
			uint64_t pos = ref_pos;
			closest_v = get_prev_vertex_with_sample(vg, idx, pos, sample_id, ref_pos,
																							sample_pos);
		}

		// start from the ref node right before the node contains sample
		closest_v = idx->find(ref_pos);
		VariantGraphVertex::sample_info sample;
		uint64_t seq_len = 0;

		if (vg->get_sample_from_vertex_if_exists(closest_v, REF, sample)){
			seq_len = ref_pos - sample.index();
			ref_pos = sample.index();
		} else {
			console->error("reference node is expected to be found!");
		}

		sample_pos = sample_pos - seq_len;

		console->debug("Closest ref node before node containing sample with sample_pos < {}: Sample pos: {}",
									 pos_x, closest_v, sample_pos);

		// DFS from such node and record seq btw [pos_x, pos_y) in sample's coordinate
		VariantGraph::VariantGraphPathIterator it = vg->find(closest_v, sample_id);
		std::string cur_ref;
		VariantGraphVertex prev_v;

		while (!it.done()) {
			Graph::vertex cur_v = (*it)->vertex_id();
			std::vector <std::string> samples_in_node;
			Variant var;
			// get next ref_pos & sample_pos
			uint64_t l = (*it)->length();
			uint64_t next_ref_pos = ref_pos + l;
			uint64_t next_sample_pos = sample_pos + l;
			std::string next_ref;

			// Check outgoing nodes, update ref_pos to the index of node containing ref
			VariantGraph::VariantGraphIterator bfs_it = vg->find((*it)->vertex_id(), 1);
			++bfs_it;
			while (!bfs_it.done()) {
				Graph::vertex v = (*bfs_it)->vertex_id();
				if (vg->get_sample_from_vertex_if_exists(v, REF, sample)) {
					next_ref_pos = sample.index();
					next_ref = vg->get_sequence(*(*bfs_it));
				}
				++bfs_it;
			}

			if (sample_pos >= pos_y) {break;}

			if (sample_pos > pos_x && vg->get_sample_from_vertex_if_exists(cur_v,
																																		 sample_id,
																																		 sample)) {
				std::string alt;
				// Insertion
				if (ref_pos == next_ref_pos) {
					cur_ref = vg->get_sequence(prev_v);
					cur_ref = "";
					console->debug("There is an insertion at (ref): {}", ref_pos);
					alt = vg->get_sequence(*(*it));
					var.var_pos = ref_pos;
				} else
					// Deletion
					if (vg->get_sample_from_vertex_if_exists(cur_v, REF, sample)) {
						console->debug("There is a deletion at (ref): {}", ref_pos);
						alt = "";
						vg->get_sample_from_vertex_if_exists(cur_v, sample_id, sample);
						var.var_pos = sample.index();
						Graph::vertex v = idx->find(ref_pos - 1); // the node before pos
						auto it = vg->find(v);
						cur_ref = vg->get_sequence(*(*it));
					} else {
						console->debug("There is a substitution at (ref): {}", ref_pos);
						alt = vg->get_sequence(*(*it));
						vg->get_sample_from_vertex_if_exists(cur_v, sample_id, sample);
						var.var_pos = sample.index();
					}

				var.alt.assign(alt);
				var.ref = cur_ref;
				get_samples((*it), vg, var.samples);

				vars.push_back(var);
			}

			cur_ref = next_ref;
			ref_pos = next_ref_pos;
			sample_pos = next_sample_pos;
			prev_v = *(*it);
			++it;
		}

		std::cout << "Number of variants get_sample_var_in_sample: " << vars.size() << '\n';
		if (print==true) {
			ofstream out;
			out.open(outfile);
			for (auto var : vars)
				print_var(&var, out);
			out.close();
		}

		return vars;
	} // get_sample_var_in_sample()


	/* ---------------------------------------------------------------------------- Given the position range and sample_id
		 Return all variants in sample_id occuring in ref coordinate [pos_x, pos_y)
		 */
	std::vector <Variant> get_sample_var_in_ref ( VariantGraph *vg, Index *idx,
																								const uint64_t pos_x,
																								const uint64_t pos_y,
																								const std::string sample_id,
																								bool print=false, std::string outfile="")
	{
		// traverse in sample's coordinate from pos_x to pos_y
		// report if the sample's node is a variant node
		std::vector <Variant> vars;
		uint64_t ref_pos;
		uint64_t sample_pos;
		Graph::vertex closest_v = get_prev_vertex_with_sample(vg, idx, pos_x,
																													sample_id,
																													ref_pos, sample_pos);
		VariantGraphVertex::sample_info sample;

		// DFS from such node and record seq btw [pos_x, pos_y) in sample's coordinate
		VariantGraph::VariantGraphPathIterator it = vg->find(closest_v, sample_id);
		std::string cur_ref;
		VariantGraphVertex prev_v;

		while (!it.done()) {
			Graph::vertex cur_v = (*it)->vertex_id();
			std::vector <std::string> samples_in_node;
			Variant var;
			// get next ref_pos & sample_pos
			uint64_t l = (*it)->length();
			uint64_t next_ref_pos = ref_pos + l;
			std::string next_ref;

			// Check outgoing nodes, update ref_pos to the index of node containing ref
			VariantGraph::VariantGraphIterator bfs_it = vg->find((*it)->vertex_id(), 1);
			++bfs_it;
			while (!bfs_it.done()) {
				Graph::vertex v = (*bfs_it)->vertex_id();
				if (vg->get_sample_from_vertex_if_exists(v, REF, sample)) {
					next_ref_pos = sample.index();
					next_ref = vg->get_sequence(*(*bfs_it));
				}
				++bfs_it;
			}

			if (ref_pos >= pos_y) {break;}

			if (ref_pos >= pos_x && vg->get_sample_from_vertex_if_exists(cur_v,
																																	 sample_id,
																																	 sample)) {
				std::string alt;
				// Insertion
				if (ref_pos == next_ref_pos) {
					cur_ref = vg->get_sequence(prev_v);
					cur_ref = "";
					console->debug("There is an insertion at (ref): {}", ref_pos);
					alt = vg->get_sequence(*(*it));
					var.var_pos = ref_pos - 1;
				} else
					// Deletion
					if (vg->get_sample_from_vertex_if_exists(cur_v, REF, sample)) {
						console->debug("There is a deletion at (ref): {}", ref_pos - 1);
						alt = "";
						Graph::vertex v = idx->find(ref_pos - 1); // the node before pos
						auto it = vg->find(v);
						cur_ref = vg->get_sequence(*(*it));
						if (vg->get_sample_from_vertex_if_exists(v, REF, sample)) {
							var.var_pos = sample.index();
						}
					} else {
						console->debug("There is a substitution at (ref): {}", ref_pos);
						alt = vg->get_sequence(*(*it));
						var.var_pos = ref_pos;
					}

				var.alt.assign(alt);
				var.ref = cur_ref;
				get_samples((*it), vg, var.samples);

				vars.push_back(var);
			}

			cur_ref = next_ref;
			ref_pos = next_ref_pos;
			prev_v = *(*it);
			++it;
		}

		std::cout << "Number of variants get_sample_var_in_ref: " << vars.size() << '\n';
		if (print==true) {
			ofstream out;
			out.open(outfile);
			for (auto var : vars)
				print_var(&var, out);
			out.close();
		}

		return vars;
	} // get_sample_var_in_ref()



	/* ---------------------------------------------------------------------------- Given the position range
		 Return all variants  occuring in ref coordinate [pos_x, pos_y)
		 */
	std::vector <Variant> get_var_in_ref ( VariantGraph *vg, Index *idx,
																				 const uint64_t pos_x,
																				 const uint64_t pos_y,
																				 bool print=false, std::string outfile="") {
		// traverse in sample's coordinate from pos_x to pos_y
		// report if the sample's node is a variant node
		// use next_variant_in_ref()
		std::vector <Variant> vars;
		uint64_t cur_pos = pos_x;
		while (cur_pos < pos_y)
		{
			uint64_t next_pos;
			if (next_variant_in_ref (vg, idx, cur_pos, vars, next_pos))
			{
				cur_pos = next_pos;
				if (cur_pos >= pos_y) {
					break;
				}
			} else {
				break;
			}
		}

		std::cout << "Number of variants get_var_in_ref: " << vars.size() << '\n';
		if (print==true) {
			ofstream out;
			out.open(outfile);
			for (auto var : vars)
				print_var(&var, out);
			out.close();
		}

		return vars;
	} // get_sample_var_in_ref



	/* ----------------------------------------------------------------------------
	// Given a variant (ref & alt) and the position
	// Return samples has such variant in the reference coordinate
	*/
	std::vector<std::pair<std::string, std::string>>
		samples_has_var(VariantGraph *vg, Index *idx, const uint64_t pos, const
										std::string ref, const std::string alt, bool print=false,
										std::string outfile="")
		{
			std::vector<Variant> vars;
			uint64_t next_pos;
			std::vector<std::pair<std::string, std::string>> samples;

			next_variant_in_ref(vg, idx, pos, vars, next_pos);
			for  (auto var : vars) {
				if (var.ref == ref && var.var_pos == pos && var.alt == alt) {
					samples.assign(var.samples.begin(),
												 var.samples.end());

					if (print == true) {
						ofstream out;
						out.open(outfile);
						for (auto i = samples.begin(); i != samples.end(); ++i)
						{
							out << i->first << ' ' << i->second;
						}
						out << std::endl;
						out.close();
					}
					return samples;
				}
			}
			console->error("There is no such variant!");

			return samples;
		} // samples_has_var()

}

#endif // __QUERY_H__
