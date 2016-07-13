#include "circularLinkedList.hpp"

namespace mantis
{
	CircularLinkedList::CircularLinkedList(int user_k, int user_minBases, int user_maxBases, int user_minRepeats, const string& user_ignoreChar)
	{
		k = user_k;
		minBases = user_minBases;
		maxBases = user_maxBases;
		minRepeats = user_minRepeats;
		curSeqLength = 0;
		numIgnoreChars = 0;
		ignoreChar = user_ignoreChar;
		lastSeqLength = 0;
		repeatCount = 0;
		first = new Node();
		cur = first;
		for (int i = 1; i < k; i++)
		{
			Node* next = new Node();
			cur->setNext(next);
			cur = next;
		}
		cur->setNext(first);
	}

	CircularLinkedList::~CircularLinkedList()
	{
		cur = first;
		for (int i = 0; i < k; i++)
		{
			Node* next = cur->getNext();
			delete cur;
			cur = next;
		}
	}

	bool CircularLinkedList::isCharInIgnore(const char& c) const
	{
		std::size_t pos = ignoreChar.find(c);
		if (pos != std::string::npos) return true;
		else return false;
	}

	bool CircularLinkedList::seqHasIgnoreChrs() const
	{
		return (numIgnoreChars > 0);
	}

	//this is a potential bottleneck - try to call this only when absolutely necessary
	void CircularLinkedList::pullKmer()
	{
		//using C strings proved to be slower than std::string +=, idk why
	/*
		Node* current = first;
		char kmerCstr[k+1];
		kmerCstr[k] = '\0';
		for (int i = 0; i < k; i++)
		{
			kmerCstr[i] = current->getElement();
			current = current->getNext();
		}
		kmer = kmerCstr;
	*/
		Node* current = first;
		kmer = "";
		for (int i = 0; i < k; i++)
		{
			kmer += current->getElement();
			current = current->getNext();
		}
	}

	bool CircularLinkedList::next(const char& c)
	{
		bool result = false;
		if (curSeqLength < k) // can only happen when seeing first k elements
		{
			bool newInIgnoreChrs = isCharN(c);
			if (newInIgnoreChrs) numIgnoreChars++;
			cur->setElement(c);
			curSeqLength++;
		}
		else
		{
			if (c == cur->getElement()) curSeqLength++;
			else
			{
				result = endSequence();
				bool curInIgnoreChrs = isCharN(cur->getElement());
				bool newInIgnoreChrs = isCharN(c);
				if (curInIgnoreChrs && !newInIgnoreChrs) numIgnoreChars--;
				else if (!curInIgnoreChrs && newInIgnoreChrs) numIgnoreChars++;
				cur->setElement(c);
				curSeqLength = k;
				first = cur->getNext();
			}
		}
		cur = cur->getNext();
		return result;
	}

	void CircularLinkedList::nextNoResultCheck(const char& c)
	{
		if (curSeqLength < k) // can only happen when seeing first k elements
		{
			bool newInIgnoreChrs = isCharN(c);
			if (newInIgnoreChrs) numIgnoreChars++;
			cur->setElement(c);
			curSeqLength++;
		}
		else
		{
			if (c == cur->getElement()) curSeqLength++;
			else
			{
				bool curInIgnoreChrs = isCharN(cur->getElement());
				bool newInIgnoreChrs = isCharN(c);
				if (curInIgnoreChrs && !newInIgnoreChrs) numIgnoreChars--;
				else if (!curInIgnoreChrs && newInIgnoreChrs) numIgnoreChars++;
				cur->setElement(c);
				curSeqLength = k;
				first = cur->getNext();
			}
		}
		cur = cur->getNext();
	}

	int CircularLinkedList::getK() const
	{
		return k;
	}

	bool CircularLinkedList::endSequence()
	{
		bool result = false;
		if (!seqHasIgnoreChrs() && curSeqLength >= minBases && curSeqLength <= maxBases)
		{
			repeatCount = curSeqLength / k;
			if (repeatCount >= minRepeats)
			{
				result = true;
				lastSeqLength = curSeqLength;
				pullKmer();
			}
		}
		return result;
	}

	unsigned int CircularLinkedList::getLastLength() const
	{
		return lastSeqLength;
	}

	string CircularLinkedList::getIgnoreChars() const
	{
		return ignoreChar;
	}

	void CircularLinkedList::setIgnoreChars(const string& user_ignoreChar)
	{
		ignoreChar = user_ignoreChar;
		numIgnoreChars = 0;
		Node* current = first;
		for (int i = 0; i < k; i++)
		{
			if (isCharInIgnore(current->getElement())) numIgnoreChars++;
			current = current->getNext();
		}
	}

	string CircularLinkedList::getKmer() const
	{
		return kmer;
	}

	int CircularLinkedList::getRepeatCount() const
	{
		return repeatCount;
	}
} //end namespace mantis
