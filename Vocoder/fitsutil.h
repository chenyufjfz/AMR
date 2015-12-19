
#ifndef FITSUTILT_H
#define FITSUTILT_H


#include <vector>

using namespace std;

namespace CCfits
{

		// vector to vector conversion.

		template <typename S, typename T>
		void fill(const vector<T>& inArray, vector<S>& outArray, size_t first, size_t number)
		{
			// vector to vector assign. stdlib takes care of deletion.
			size_t last = number + first -1;
			outArray.resize(number);
			for (size_t j = first; j <= last; ++j)
			{
				outArray[j - first] = static_cast<S>(inArray[j]);
			}
		}


		template <typename S, typename T>
		void fill(const vector<T>& inArray, vector<S>& outArray)
		{
			fill(inArray, outArray, 0, inArray.size());
		}

		template <typename S, typename T>
		void fill(const T * inArray, vector<S>& outArray, size_t first, size_t number)
		{
			// vector to vector assign. stdlib takes care of deletion.
			size_t last = number + first -1;
			outArray.resize(number);
			for (size_t j = first; j <= last; ++j)
			{
				outArray[j - first] = static_cast<S>(inArray[j]);
			}
		}

		template <typename S, typename T>
		void fill(const vector<T>& inArray, S * outArray, size_t first, size_t number)
		{
			// vector to vector assign. stdlib takes care of deletion.
			size_t last = number + first -1;
			for (size_t j = first; j <= last; ++j)
			{
				outArray[j - first] = static_cast<S>(inArray[j]);
			}
		}

		template <typename S, typename T>
		void fill(const vector<T>& inArray, S * outArray)
		{
			fill(inArray, outArray, 0, inArray.size());
		}

} // namespace CCfits


#endif
