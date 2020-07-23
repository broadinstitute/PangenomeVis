/**
 * Quantity Counter wrapper.
 *
 * @author Htmlstream
 * @version 1.0
 *
 */
;(function ($) {
	'use strict';

	$.HSCore.components.HSQantityCounter = {
		/**
		 *
		 *
		 * @var Object _baseConfig
		 */
		_baseConfig: {
			minVal: 0,
			maxVal: Infinity,
			beforeIncrease: function () {},
			afterIncrease: function () {},
			stopIncrease: function() {},

			beforeDecrease: function () {},
			afterDecrease: function () {},
			stopDecrease: function() {},

			beforeChange: function () {},
			afterChange: function () {}
		},

		/**
		 *
		 *
		 * @var jQuery pageCollection
		 */
		pageCollection: $(),

		/**
		 * Initialization of Count quantity wrapper.
		 *
		 * @param String selector (optional)
		 * @param Object config (optional)
		 *
		 * @return jQuery pageCollection - collection of initialized items.
		 */

		init: function (selector, config) {

			this.collection = selector && $(selector).length ? $(selector) : $();
			if (!$(selector).length) return;

			this.config = config && $.isPlainObject(config) ?
				$.extend({}, this._baseConfig, config) : this._baseConfig;

			this.config.itemSelector = selector;

			this.initCountQty();

			return this.pageCollection;

		},

		initCountQty: function () {
			//Variables
			var $self = this,
				config = $self.config,
				collection = $self.pageCollection;

			//Actions
			this.collection.each(function (i, el) {
				//Variables
				var $this = $(el),
					$plus = $this.find('.js-plus'),
					$minus = $this.find('.js-minus'),
					$result = $this.find('.js-result'),
					minVal = $this.data('min-val') ? $this.data('min-val') : config.minVal,
					maxVal = $this.data('max-val') ? $this.data('max-val') - 2 : config.maxVal,
					resultVal = parseInt($result[0].attributes.value !== undefined ? $result.val() : $result.text());

				$plus.on('click', function (e) {
					e.preventDefault();

					if (config.stopIncrease()) return;

					config.beforeIncrease();
					config.beforeChange();

					if (resultVal <= (maxVal + 1)) {

						resultVal += 1;

						if ($result[0].attributes.value) {

							$result.val(resultVal);

						} else {

							$result.text(resultVal);

						}

					}

					config.afterIncrease();
					config.afterChange();
				});

				$minus.on('click', function (e) {
					e.preventDefault();

					if (config.stopDecrease()) return;

					config.beforeDecrease();
					config.beforeChange();

					if (resultVal >= (minVal + 1)) {
						resultVal -= 1;

						if ($result[0].attributes.value) {

							$result.val(resultVal);

						} else {

							$result.text(resultVal);

						}
					} else {

						return false;

					}

					config.afterDecrease();
					config.afterChange();
				});

				//Actions
				collection = collection.add($this);
			});
		}

	};

})(jQuery);