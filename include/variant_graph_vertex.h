/*
 * ============================================================================
 *
 *       Filename:  variant_graph_vertex.h
 *
 *         Author:  Prashant Pandey (), ppandey2@cs.cmu.edu
 *   Organization:  Carnegie Mellon University
 *
 * ============================================================================
 */

#ifndef __VARIANT_GRAPH_VERTEX_H__
#define __VARIANT_GRAPH_VERTEX_H__

#include <string>
#include <vector>

namespace variantdb {

	// VariantGraphVertex class to store information about each vertex.
	// Each vertex contains a @vertex_id, @offset, and @length.
	// Multiple samples can be attached to a vertex.
	// A sample_info contains an @index, @sample_id, and genotype info.

	class VariantGraphVertex {
		public:
			class sample_info {
				public:
					sample_info() : index_(0), sample_id_(""), gt_1_(false),
					gt_2_(false) {}
					sample_info(uint32_t idx, const std::string id, bool gt1, bool gt2)
						: index_(idx), sample_id_(id), gt_1_(gt1), gt_2_(gt2) {}

					uint32_t index(void) const { return index_; }
					const std::string& sample_id(void) const { return sample_id_; }
					bool gt_1(void) const { return gt_1_; }
					bool gt_2(void) const { return gt_2_; }

				private:
					uint32_t index_;
					std::string sample_id_;
					// two bits for genotype info.
					bool gt_1_;
					bool gt_2_;
			};

			VariantGraphVertex(uint32_t id, uint32_t off, uint32_t len, const
												 sample_info& s) : vertex_id_(id), offset_(off),
			length_(len) {
				s_info_.emplace_back(s);
			}

			~VariantGraphVertex() {
				s_info_.clear();
			}

			void add_sample(const sample_info& s) {
				s_info_.emplace_back(s);
			}

			void set_length(uint32_t len) {
				length_ = len;
			}

			uint32_t vertex_id(void) const { return vertex_id_; }
			uint32_t offset(void) const { return offset_; }
			uint32_t length(void) const { return length_; }
			uint32_t s_info_size(void) const { return s_info_.size(); }
			const sample_info& s_info(uint32_t index) const {
				return s_info_[index];
			} 

		private:
			uint32_t vertex_id_;
			uint32_t offset_;
			uint32_t length_;
			std::vector<sample_info> s_info_;
	};

	class VariantGraphVertexList {
		public:
			VariantGraphVertexList() {}
			VariantGraphVertexList(const std::string file);

			~VariantGraphVertexList() {
				vertexes.clear();
			}

			const VariantGraphVertex* add_vertex(const VariantGraphVertex& v) {
				vertexes.emplace_back(v);
				return &vertexes[vertexes.size() - 1];
			}

			const VariantGraphVertex& vertex(uint32_t index) const {
				return vertexes[index];
			}

			VariantGraphVertex* mutable_vertex(uint32_t index) {
				return &vertexes[index];
			};

			void serialize(const std::string file);

		private:
			std::vector<VariantGraphVertex> vertexes;
	};

}

#endif // __VARIANT_GRAPH_VERTEX_H__
